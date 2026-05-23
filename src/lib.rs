//! Streaming PCR/optical-duplicate marking, a Rust port of `samtools markdup`
//! (MIT) default template mode.
//!
//! Input is coordinate-sorted and pre-processed by `samtools fixmate -m`, so
//! every paired record carries `MC` (mate CIGAR) and `ms` (mate score). markdup
//! streams in coordinate order keeping only a bounded window of recent reads:
//! a read held in the buffer is finalized and flushed once the current read's
//! coordinate has advanced past it by more than `max_length` bases or the
//! reference changes (samtools `bam_mark_duplicates`, the buffer-trim loop).
//! When a record leaves the window its single/pair hash slots are removed, so
//! memory stays bounded by the window span rather than the file size — the
//! whole point of the rewrite over the prior full-buffer version.
//!
//! Only the 0x400 duplicate flag bit is edited; seq/qual/cigar/name pass
//! through byte-for-byte via the [`rsomics_bamio::raw`] path.

mod key;

use std::collections::HashMap;
use std::collections::VecDeque;
use std::io::Write;
use std::num::NonZero;
use std::path::Path;

use noodles::bam;
use rsomics_bamio::raw::{self, FLAG_DUPLICATE, RawRecord};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use key::{PairKey, SingleKey, make_pair_key, make_single_key};

const FLAG_PAIRED: u16 = 0x1;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QCFAIL: u16 = 0x200;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// samtools default minimum base quality counted toward a read's score.
const MD_MIN_QUALITY: u8 = 15;

/// samtools default sliding-window span in bases (`-l`, default 300).
const DEFAULT_MAX_LENGTH: i64 = 300;

#[derive(Debug, Default, Clone, Serialize)]
pub struct MarkdupStats {
    pub total: u64,
    pub duplicates_marked: u64,
    pub duplicates_removed: u64,
}

#[derive(Debug, Clone, Default)]
pub struct MarkdupOpts {
    pub remove: bool,
}

/// Sum of base qualities >= 15 (samtools `calc_score`).
fn calc_score(r: &RawRecord) -> i64 {
    r.quality_scores()
        .iter()
        .filter(|&&q| q >= MD_MIN_QUALITY)
        .map(|&q| i64::from(q))
        .sum()
}

/// Decode the `ms` mate-score tag, accepting any BAM integer subtype
/// (`c`/`C`/`s`/`S`/`i`/`I`) exactly as samtools `bam_aux2i` does. Absent or
/// non-integer `ms` is a hard error: the file was not fixmate-m'd.
fn mate_score(r: &RawRecord) -> Result<i64> {
    let ty = r.aux_type(*b"ms").ok_or_else(no_ms_error)?;
    let v = r.aux_value(*b"ms").ok_or_else(no_ms_error)?;
    let val = match ty {
        b'c' => i64::from(v[0] as i8),
        b'C' => i64::from(v[0]),
        b's' => i64::from(i16::from_le_bytes([v[0], v[1]])),
        b'S' => i64::from(u16::from_le_bytes([v[0], v[1]])),
        b'i' => i64::from(i32::from_le_bytes([v[0], v[1], v[2], v[3]])),
        b'I' => i64::from(u32::from_le_bytes([v[0], v[1], v[2], v[3]])),
        _ => {
            return Err(RsomicsError::InvalidInput(
                "ms tag is not an integer type. Run samtools fixmate -m first.".into(),
            ));
        }
    };
    Ok(val)
}

fn no_ms_error() -> RsomicsError {
    RsomicsError::InvalidInput("no ms score tag. Run samtools fixmate -m on the file first.".into())
}

fn no_mc_error() -> RsomicsError {
    RsomicsError::InvalidInput("no MC tag. Run samtools fixmate -m on the file first.".into())
}

/// A read with a mate: paired, mate mapped, mate reference/pos present.
/// Mirrors samtools `has_mate`.
fn has_mate(r: &RawRecord) -> bool {
    r.flags() & FLAG_PAIRED != 0
        && r.flags() & FLAG_MATE_UNMAPPED == 0
        && !(r.mate_reference_sequence_id() == -1 && r.mate_alignment_start() == -1)
}

/// A buffered read awaiting flush. Holds its raw bytes, the window position
/// (the single key's `this_coord`), its raw tid, and the hash keys it occupies
/// so they can be removed on flush.
struct Buffered {
    record: RawRecord,
    is_dup: bool,
    pos: i64,
    tid: i32,
    /// single-hash key, present iff this read registered a slot it still owns.
    single_key: Option<SingleKey>,
    /// pair-hash key, present iff this read registered a pair slot it still owns.
    pair_key: Option<PairKey>,
}

/// Value stored in the single/pair hashes: the buffer index of the occupying
/// read. Generation-tagged buffer indices are unnecessary because flushed reads
/// remove their own slots before the index could be reused.
type SingleHash = HashMap<SingleKey, usize>;
type PairHash = HashMap<PairKey, usize>;

/// Streaming duplicate marker. Owns the sliding-window buffer and the two
/// position hashes; emits finalized records to `writer` in input order.
struct Marker<'a, W: Write> {
    writer: &'a mut bam::io::Writer<W>,
    opts: &'a MarkdupOpts,
    buffer: VecDeque<Buffered>,
    /// Buffer index of buffer[0], so hash values (absolute indices) map to slots.
    base: usize,
    single_hash: SingleHash,
    pair_hash: PairHash,
    max_length: i64,
    stats: MarkdupStats,
}

impl<'a, W: Write> Marker<'a, W> {
    fn new(writer: &'a mut bam::io::Writer<W>, opts: &'a MarkdupOpts) -> Self {
        Marker {
            writer,
            opts,
            buffer: VecDeque::new(),
            base: 0,
            single_hash: HashMap::new(),
            pair_hash: HashMap::new(),
            max_length: DEFAULT_MAX_LENGTH,
            stats: MarkdupStats::default(),
        }
    }

    #[inline]
    fn slot(&self, abs_idx: usize) -> &Buffered {
        &self.buffer[abs_idx - self.base]
    }

    #[inline]
    fn slot_mut(&mut self, abs_idx: usize) -> &mut Buffered {
        &mut self.buffer[abs_idx - self.base]
    }

    /// Process one record: build keys, resolve duplicates against the live
    /// window hashes, then push it onto the buffer. Excluded reads (secondary /
    /// supplementary / unmapped / QC-fail) are buffered untouched so they emit
    /// in order.
    fn process(&mut self, record: RawRecord) -> Result<()> {
        self.stats.total += 1;
        let flags = record.flags();
        let raw_pos = i64::from(record.alignment_start());
        let raw_tid = record.reference_sequence_id();
        let excluded =
            flags & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_UNMAPPED | FLAG_QCFAIL) != 0;

        if excluded {
            self.push(Buffered {
                record,
                is_dup: false,
                pos: raw_pos,
                tid: raw_tid,
                single_key: None,
                pair_key: None,
            });
            self.flush_window(raw_pos, raw_tid)?;
            return Ok(());
        }

        let single_key = make_single_key(&record);
        let window_pos = single_key.this_coord;

        let mut buffered = Buffered {
            record,
            is_dup: false,
            pos: window_pos,
            tid: raw_tid,
            single_key: None,
            pair_key: None,
        };

        let new_idx = self.base + self.buffer.len();
        let paired = has_mate(&buffered.record);

        if paired {
            let mc = buffered
                .record
                .aux_value(*b"MC")
                .ok_or_else(no_mc_error)?
                // MC is a NUL-terminated Z string; drop the terminator for parsing.
                .split(|&b| b == 0)
                .next()
                .unwrap_or(&[]);
            let pair_key = make_pair_key(&buffered.record, mc);

            self.resolve_single(&mut buffered, single_key, new_idx, true)?;
            self.resolve_pair(&mut buffered, pair_key, new_idx)?;
        } else {
            self.resolve_single(&mut buffered, single_key, new_idx, false)?;
        }

        self.push(buffered);
        self.flush_window(raw_pos, raw_tid)?;
        Ok(())
    }

    /// single_hash interaction. `paired` marks whether the incoming read has a
    /// mate (its slot wins over singletons that later collide).
    fn resolve_single(
        &mut self,
        buffered: &mut Buffered,
        single_key: SingleKey,
        new_idx: usize,
        paired: bool,
    ) -> Result<()> {
        match self.single_hash.get(&single_key).copied() {
            None => {
                self.single_hash.insert(single_key, new_idx);
                buffered.single_key = Some(single_key);
            }
            Some(occupant_idx) => {
                let occupant_paired = has_mate(&self.slot(occupant_idx).record);
                if paired {
                    // Incoming paired read collides with an existing slot.
                    if !occupant_paired {
                        // Singleton always loses to a pair; the singleton is dup
                        // and the pair takes the slot.
                        self.mark_dup(occupant_idx);
                        self.slot_mut(occupant_idx).single_key = None;
                        self.single_hash.insert(single_key, new_idx);
                        buffered.single_key = Some(single_key);
                    }
                    // occupant is itself paired: pair_hash resolves them.
                } else if occupant_paired {
                    // Incoming singleton matched against one of a pair: it's dup.
                    buffered.is_dup = true;
                    self.stats.duplicates_marked += 1;
                } else {
                    // Singleton vs singleton: higher base-qual score wins; ties
                    // keep the first (coordinate order, samtools behaviour).
                    let old_score = calc_score(&self.slot(occupant_idx).record);
                    let new_score = calc_score(&buffered.record);
                    if new_score > old_score {
                        self.mark_dup(occupant_idx);
                        self.slot_mut(occupant_idx).single_key = None;
                        self.single_hash.insert(single_key, new_idx);
                        buffered.single_key = Some(single_key);
                    } else {
                        buffered.is_dup = true;
                        self.stats.duplicates_marked += 1;
                    }
                }
            }
        }
        Ok(())
    }

    /// pair_hash interaction with the summed-mate-score tie-break.
    fn resolve_pair(
        &mut self,
        buffered: &mut Buffered,
        pair_key: PairKey,
        new_idx: usize,
    ) -> Result<()> {
        match self.pair_hash.get(&pair_key).copied() {
            None => {
                self.pair_hash.insert(pair_key, new_idx);
                buffered.pair_key = Some(pair_key);
            }
            Some(occupant_idx) => {
                let occupant_qcfail = self.slot(occupant_idx).record.flags() & FLAG_QCFAIL != 0;
                let new_qcfail = buffered.record.flags() & FLAG_QCFAIL != 0;

                let (old_score, new_score) = if occupant_qcfail != new_qcfail {
                    // The non-QC-fail read wins regardless of base scores.
                    if occupant_qcfail { (0, 1) } else { (1, 0) }
                } else {
                    let old = calc_score(&self.slot(occupant_idx).record)
                        + mate_score(&self.slot(occupant_idx).record)?;
                    let new = calc_score(&buffered.record) + mate_score(&buffered.record)?;
                    (old, new)
                };

                // Tie-break by qname: lexicographically smaller name is original
                // (samtools: strcmp(in_read, occupant) < 0 -> tie_add = +1).
                let tie_add = if new_score == old_score {
                    if buffered.record.name() < self.slot(occupant_idx).record.name() {
                        1
                    } else {
                        -1
                    }
                } else {
                    0
                };

                if new_score + tie_add > old_score {
                    // Incoming read wins: occupant becomes dup, incoming takes slot.
                    self.mark_dup(occupant_idx);
                    self.slot_mut(occupant_idx).pair_key = None;
                    self.pair_hash.insert(pair_key, new_idx);
                    buffered.pair_key = Some(pair_key);
                } else {
                    buffered.is_dup = true;
                    self.stats.duplicates_marked += 1;
                }
            }
        }
        Ok(())
    }

    /// Mark a buffered read as a duplicate, accounting for it once. A read marked
    /// from the single hash and then re-evaluated cannot double-count because
    /// the slot's `single_key`/`pair_key` are cleared on the same step.
    fn mark_dup(&mut self, abs_idx: usize) {
        let slot = self.slot_mut(abs_idx);
        if !slot.is_dup {
            slot.is_dup = true;
            self.stats.duplicates_marked += 1;
        }
    }

    fn push(&mut self, b: Buffered) {
        self.buffer.push_back(b);
    }

    /// Flush every buffered read that has fallen out of the window relative to
    /// the current read's raw `(pos, tid)`. samtools keeps a read while
    /// `read.pos + max_length > cur_pos && read.tid == cur_tid`; everything else
    /// is finalized, its hash slots removed, and written out.
    fn flush_window(&mut self, cur_pos: i64, cur_tid: i32) -> Result<()> {
        while let Some(front) = self.buffer.front() {
            let keep = front.pos + self.max_length > cur_pos
                && front.tid == cur_tid
                && (cur_tid != -1 || cur_pos != -1);
            if keep {
                break;
            }
            self.emit_front()?;
        }
        Ok(())
    }

    /// Remove the front buffered read's hash slots and write it out (or skip on
    /// remove). Advances `base`.
    fn emit_front(&mut self) -> Result<()> {
        let b = self.buffer.pop_front().expect("emit_front on empty buffer");
        self.base += 1;
        if let Some(sk) = b.single_key {
            self.single_hash.remove(&sk);
        }
        if let Some(pk) = b.pair_key {
            self.pair_hash.remove(&pk);
        }
        self.emit(b)
    }

    fn emit(&mut self, b: Buffered) -> Result<()> {
        if b.is_dup {
            if self.opts.remove {
                self.stats.duplicates_removed += 1;
                return Ok(());
            }
            let mut edited = b.record;
            edited.set_flag_bits(FLAG_DUPLICATE);
            raw::write_record(self.writer.get_mut(), &edited)?;
        } else {
            raw::write_record(self.writer.get_mut(), &b.record)?;
        }
        Ok(())
    }

    /// Flush all remaining buffered reads at end of stream.
    fn finish(mut self) -> Result<MarkdupStats> {
        while !self.buffer.is_empty() {
            self.emit_front()?;
        }
        Ok(self.stats)
    }
}

/// Stream `input` (coordinate-sorted, fixmate-m'd) and emit duplicate-marked
/// records. `output_path` of `None` writes BAM to stdout.
pub fn markdup(
    input: &Path,
    output_path: Option<&Path>,
    opts: &MarkdupOpts,
    workers: NonZero<usize>,
) -> Result<MarkdupStats> {
    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    match output_path {
        Some(path) => {
            let mut writer = rsomics_bamio::create_with_workers(path, workers)?;
            run(&mut reader, &header, &mut writer, opts)
        }
        None => {
            let mut writer = bam::io::Writer::new(std::io::stdout().lock());
            run(&mut reader, &header, &mut writer, opts)
        }
    }
}

fn run<R, W>(
    reader: &mut bam::io::Reader<R>,
    header: &noodles::sam::Header,
    writer: &mut bam::io::Writer<W>,
    opts: &MarkdupOpts,
) -> Result<MarkdupStats>
where
    R: std::io::Read,
    W: Write,
{
    writer.write_header(header).map_err(RsomicsError::Io)?;

    let mut marker = Marker::new(writer, opts);
    let mut rec = RawRecord::default();
    while raw::read_record(reader.get_mut(), &mut rec)? != 0 {
        marker.process(std::mem::take(&mut rec))?;
    }
    marker.finish()
}
