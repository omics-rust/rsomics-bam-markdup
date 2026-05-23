use std::collections::HashMap;
use std::io::Write;
use std::num::NonZero;
use std::path::Path;

use noodles::bam;
use rsomics_bamio::raw::{self, FLAG_DUPLICATE, RawRecord};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

// SAM/BAM FLAG bits (SAMv1 §1.4); duplicate (0x400) lives in rsomics-bamio.
const FLAG_SEGMENTED: u16 = 0x1;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_MATE_UNMAPPED: u16 = 0x8;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_MATE_REVERSE: u16 = 0x20;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

// CIGAR op codes (BAM packed encoding, low nibble): M=0 I=1 D=2 N=3 S=4 H=5 P=6 ==7 X=8.
const CIGAR_MATCH: u8 = 0;
const CIGAR_DELETION: u8 = 2;
const CIGAR_SKIP: u8 = 3;
const CIGAR_SOFT_CLIP: u8 = 4;
const CIGAR_SEQ_MATCH: u8 = 7;
const CIGAR_SEQ_MISMATCH: u8 = 8;

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

fn is_reverse(r: &RawRecord) -> bool {
    r.flags() & FLAG_REVERSE != 0
}

fn is_mate_reverse(r: &RawRecord) -> bool {
    r.flags() & FLAG_MATE_REVERSE != 0
}

/// 1-based alignment start (BAM stores 0-based; -1 → 0, matching the unmapped
/// sentinel handling of the prior noodles-based path).
fn alignment_start_1based(r: &RawRecord) -> i64 {
    let pos = r.alignment_start();
    if pos < 0 { 0 } else { i64::from(pos) + 1 }
}

/// (reference-consuming span, leading soft-clip, trailing soft-clip) from CIGAR.
fn span_and_clips(r: &RawRecord) -> (i64, i64, i64) {
    let ops: Vec<(u8, u32)> = r.cigar_ops().collect();
    let leading = ops
        .first()
        .filter(|(k, _)| *k == CIGAR_SOFT_CLIP)
        .map_or(0, |(_, len)| i64::from(*len));
    let trailing = ops
        .last()
        .filter(|(k, _)| *k == CIGAR_SOFT_CLIP)
        .map_or(0, |(_, len)| i64::from(*len));
    let span: i64 = ops
        .iter()
        .filter(|(k, _)| {
            matches!(
                *k,
                CIGAR_MATCH | CIGAR_DELETION | CIGAR_SKIP | CIGAR_SEQ_MATCH | CIGAR_SEQ_MISMATCH
            )
        })
        .map(|(_, len)| i64::from(*len))
        .sum();
    (span, leading, trailing)
}

/// Unclipped 5' coordinate (start): forward = pos - leading_clip.
fn unclipped_start(r: &RawRecord) -> i64 {
    let (_, lead, _) = span_and_clips(r);
    alignment_start_1based(r) - lead
}

/// Unclipped end coordinate: pos + ref_span + trailing_clip.
fn unclipped_end(r: &RawRecord) -> i64 {
    let (span, _, trail) = span_and_clips(r);
    alignment_start_1based(r) + span + trail
}

/// Unclipped 5' coordinate: forward = unclipped_start; reverse = unclipped_end.
fn unclipped_5p(r: &RawRecord) -> i64 {
    if is_reverse(r) {
        unclipped_end(r)
    } else {
        unclipped_start(r)
    }
}

fn tid(r: &RawRecord) -> i64 {
    i64::from(r.reference_sequence_id())
}

fn mate_tid(r: &RawRecord) -> i64 {
    i64::from(r.mate_reference_sequence_id())
}

fn base_qual_sum(r: &RawRecord) -> i64 {
    // samtools markdup sums only base qualities >= 15 (its default threshold).
    r.quality_scores()
        .iter()
        .filter(|&&q| q >= 15)
        .map(|&q| i64::from(q))
        .sum()
}

/// Single-end key: (tid, unclipped_5', is_reverse).
/// Mirrors samtools make_single_key (without barcode/RG support).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct SingleKey {
    tid: i64,
    coord: i64,
    rev: bool,
}

/// Pair key from the perspective of one read, matching samtools template-mode make_pair_key.
///
/// samtools stores `(this_coord, other_coord, leftmost, orientation)` from the
/// processing read's perspective — so the two mates of a template produce
/// DIFFERENT keys (left-read perspective vs right-read perspective).  Only
/// reads from distinct templates at identical positions share a key, which is
/// the actual duplicate condition.
///
/// Orientation encoding (mirrors samtools O_FF/O_RR/O_FR/O_RF):
///   0 = O_FF, 1 = O_RR, 2 = O_FR, 3 = O_RF
///
/// leftmost: true = this read is the leftmost end (R_LE), false = R_RI.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct PairKey {
    this_tid: i64,
    this_coord: i64,
    other_tid: i64,
    other_coord: i64,
    leftmost: bool,
    orientation: u8,
}

fn single_key_for(r: &RawRecord) -> SingleKey {
    SingleKey {
        tid: tid(r),
        coord: unclipped_5p(r),
        rev: is_reverse(r),
    }
}

/// Build the pair key from this read's perspective, given the mate's unclipped_5' coord.
/// Mirrors samtools template-mode make_pair_key (no barcode/RG).
fn pair_key_for(r: &RawRecord, mate_5p: i64) -> PairKey {
    let this_ref = tid(r);
    let other_ref = mate_tid(r);
    let this_rev = is_reverse(r);
    let mate_rev = is_mate_reverse(r);

    let this_start = unclipped_start(r);
    let this_end = unclipped_end(r);

    // Determine leftmost in template mode (mirroring samtools):
    let leftmost = if this_ref != other_ref {
        this_ref < other_ref
    } else if this_rev == mate_rev {
        // Same orientation (FF or RR):
        if !this_rev {
            this_start <= mate_5p // forward: compare starts
        } else {
            this_end <= mate_5p // reverse: compare ends
        }
    } else {
        // Mixed orientation (FR or RF):
        if this_rev {
            this_end <= mate_5p // reverse this: compare this_end vs mate_start
        } else {
            this_start <= mate_5p // forward this: compare this_start vs mate_end
        }
    };

    let orientation: u8 = match (this_rev, mate_rev) {
        (false, false) => 0, // O_FF
        (true, true) => 1,   // O_RR
        (false, true) => 2,  // O_FR
        (true, false) => 3,  // O_RF
    };

    // Coordinate to store: forward read uses unclipped_start; reverse uses unclipped_end.
    let this_coord = if this_rev { this_end } else { this_start };
    // Other coord: mate_5p is already the mate's unclipped_5' (start for fwd, end for rev).
    let other_coord = mate_5p;

    PairKey {
        this_tid: this_ref,
        this_coord,
        other_tid: other_ref,
        other_coord,
        leftmost,
        orientation,
    }
}

/// Entry in the single-hash: record index and whether the occupant is a PE read.
#[derive(Debug, Clone, Copy)]
struct SingleEntry {
    idx: usize,
    is_paired: bool,
}

/// Mark per-record `is_dup` flags using the position-keyed two-hash detection
/// that mirrors samtools markdup. Returns a bitmap parallel to `records`.
fn detect_duplicates(records: &[RawRecord]) -> Vec<bool> {
    // --- Pass 1: resolve mate unclipped-5' coordinates for paired reads. ---
    //
    // samtools markdup requires the MC (mate CIGAR) tag from fixmate to compute
    // the mate's unclipped position.  We instead do a name-grouping pass over
    // the already-loaded records, giving us both ends without requiring fixmate.
    // Mate's resolved unclipped-5' coord per record index, or `None` for reads
    // with no mate to pair. Indexed by record position to avoid a second hash.
    let mut mate_5p: Vec<Option<i64>> = vec![None; records.len()];
    {
        // Name keys borrow from `records` (immutable here) — no per-read clone.
        let mut name_to_idxs: HashMap<&[u8], Vec<usize>> = HashMap::new();
        for (i, r) in records.iter().enumerate() {
            let f = r.flags();
            if f & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
                continue;
            }
            if f & FLAG_SEGMENTED != 0 && f & FLAG_MATE_UNMAPPED == 0 {
                name_to_idxs.entry(r.name()).or_default().push(i);
            }
        }
        for idxs in name_to_idxs.values() {
            if idxs.len() != 2 {
                // Chimeric / multi-segment: treat each end as single-end.
                continue;
            }
            let (i0, i1) = (idxs[0], idxs[1]);
            mate_5p[i0] = Some(unclipped_5p(&records[i1]));
            mate_5p[i1] = Some(unclipped_5p(&records[i0]));
        }
    }

    // Base-quality score per record, computed once (the SE-vs-SE collision rule
    // needs the incumbent's score, so caching avoids a second CIGAR-free walk).
    let scores: Vec<i64> = records.iter().map(base_qual_sum).collect();

    // --- Pass 2: position-keyed duplicate detection. ---
    //
    // Two hash tables mirror samtools' single_hash and pair_hash:
    //
    //   single_hash: SingleKey → SingleEntry
    //     Holds every primary mapped read keyed by its own unclipped 5' position
    //     + orientation.  A PE read that occupies a slot beats any SE read that
    //     later collides (PE always wins over SE at the same position).
    //
    //   pair_hash: PairKey → (score, idx)
    //     Holds paired reads keyed by the read-perspective pair key.  Because the
    //     key is NOT symmetric (leftmost flag + this vs other coord), two mates of
    //     the same template produce different keys and never collide — only reads
    //     from distinct templates at identical positions collide.
    //
    // Collision rules (matching samtools default behaviour):
    //   SE vs SE in single_hash  → higher base-qual score wins; ties keep first.
    //   PE vs SE in single_hash  → PE always wins; SE marked dup.
    //   SE vs PE in single_hash  → PE already in slot wins; SE marked dup.
    //   PE vs PE in pair_hash    → higher per-read score wins; ties keep first.

    let mut single_hash: HashMap<SingleKey, SingleEntry> = HashMap::new();
    let mut pair_hash: HashMap<PairKey, (i64, usize)> = HashMap::new();
    let mut is_dup = vec![false; records.len()];

    for (i, record) in records.iter().enumerate() {
        let f = record.flags();
        if f & (FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY) != 0 {
            continue;
        }

        let sk = single_key_for(record);
        let score = scores[i];
        let is_paired =
            f & FLAG_SEGMENTED != 0 && f & FLAG_MATE_UNMAPPED == 0 && mate_5p[i].is_some();

        // --- single_hash: every read registers its own-end key ---
        match single_hash.get(&sk).copied() {
            None => {
                single_hash.insert(sk.clone(), SingleEntry { idx: i, is_paired });
            }
            Some(existing) if existing.is_paired && !is_paired => {
                // SE collides with an existing PE entry: PE wins, SE is dup.
                is_dup[i] = true;
            }
            Some(existing) if !existing.is_paired && is_paired => {
                // PE collides with an existing SE entry: PE wins, SE is dup.
                is_dup[existing.idx] = true;
                single_hash.insert(
                    sk.clone(),
                    SingleEntry {
                        idx: i,
                        is_paired: true,
                    },
                );
            }
            Some(existing) if !existing.is_paired && !is_paired => {
                // SE vs SE: higher score wins; ties keep the first (coordinate order).
                let old_score = scores[existing.idx];
                if score > old_score {
                    is_dup[existing.idx] = true;
                    single_hash.insert(
                        sk.clone(),
                        SingleEntry {
                            idx: i,
                            is_paired: false,
                        },
                    );
                } else {
                    is_dup[i] = true;
                }
            }
            Some(_) => {
                // PE vs PE in single_hash: pair_hash resolves this below.
            }
        }

        // --- pair_hash: only for reads with a resolved mate ---
        if let Some(mate_coord) = mate_5p[i] {
            let pk = pair_key_for(record, mate_coord);

            match pair_hash.get(&pk).copied() {
                None => {
                    pair_hash.insert(pk, (score, i));
                }
                Some((old_score, old_idx)) => {
                    // Two reads at identical pair positions: higher score wins.
                    if score > old_score {
                        is_dup[old_idx] = true;
                        pair_hash.insert(pk, (score, i));
                    } else {
                        is_dup[i] = true;
                    }
                }
            }
        }
    }

    is_dup
}

/// Read every record raw, detect duplicates, then for each duplicate set the
/// 0x400 flag in place and emit the raw bytes — seq/qual/cigar/name are never
/// decoded or re-encoded. `output_path` of `None` writes BAM to stdout.
pub fn markdup(
    input: &Path,
    output_path: Option<&Path>,
    opts: &MarkdupOpts,
    workers: NonZero<usize>,
) -> Result<MarkdupStats> {
    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    let mut records: Vec<RawRecord> = Vec::new();
    let mut rec = RawRecord::default();
    while raw::read_record(reader.get_mut(), &mut rec)? != 0 {
        records.push(std::mem::take(&mut rec));
    }

    let is_dup = detect_duplicates(&records);

    match output_path {
        Some(path) => {
            let mut writer = rsomics_bamio::create_with_workers(path, workers)?;
            write_records(&mut writer, &header, &records, &is_dup, opts)
        }
        None => {
            let mut writer = bam::io::Writer::new(std::io::stdout().lock());
            write_records(&mut writer, &header, &records, &is_dup, opts)
        }
    }
}

fn write_records<W: Write>(
    writer: &mut bam::io::Writer<W>,
    header: &noodles::sam::Header,
    records: &[RawRecord],
    is_dup: &[bool],
    opts: &MarkdupOpts,
) -> Result<MarkdupStats> {
    writer.write_header(header).map_err(RsomicsError::Io)?;

    let mut stats = MarkdupStats {
        total: records.len() as u64,
        ..Default::default()
    };

    for (i, record) in records.iter().enumerate() {
        if is_dup[i] {
            if opts.remove {
                stats.duplicates_removed += 1;
                continue;
            }
            stats.duplicates_marked += 1;
            let mut edited = record.clone();
            edited.set_flag_bits(FLAG_DUPLICATE);
            raw::write_record(writer.get_mut(), &edited)?;
        } else {
            raw::write_record(writer.get_mut(), record)?;
        }
    }

    Ok(stats)
}
