use std::collections::HashMap;
use std::io::Write;
use std::num::NonZero;
use std::path::Path;

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as AlnWrite;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::RecordBuf;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

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

/// (reference-consuming span, leading soft-clip, trailing soft-clip) from CIGAR.
fn span_and_clips(record: &bam::Record) -> (i64, i64, i64) {
    let ops: Vec<_> = record
        .cigar()
        .iter()
        .filter_map(std::result::Result::ok)
        .collect();
    let leading = ops
        .first()
        .filter(|o| o.kind() == Kind::SoftClip)
        .map_or(0, |o| o.len() as i64);
    let trailing = ops
        .last()
        .filter(|o| o.kind() == Kind::SoftClip)
        .map_or(0, |o| o.len() as i64);
    let span: i64 = ops
        .iter()
        .filter(|o| {
            matches!(
                o.kind(),
                Kind::Match
                    | Kind::Deletion
                    | Kind::Skip
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
            )
        })
        .map(|o| o.len() as i64)
        .sum();
    (span, leading, trailing)
}

/// Unclipped 5' coordinate (start): forward = pos - leading_clip; same for unclipped_start.
fn unclipped_start(record: &bam::Record) -> i64 {
    let pos = record
        .alignment_start()
        .and_then(std::result::Result::ok)
        .map_or(0, |p| p.get() as i64);
    let (_, lead, _) = span_and_clips(record);
    pos - lead
}

/// Unclipped end coordinate: pos + ref_span + trailing_clip.
fn unclipped_end(record: &bam::Record) -> i64 {
    let pos = record
        .alignment_start()
        .and_then(std::result::Result::ok)
        .map_or(0, |p| p.get() as i64);
    let (span, _, trail) = span_and_clips(record);
    pos + span + trail
}

/// Unclipped 5' coordinate: forward = unclipped_start; reverse = unclipped_end.
fn unclipped_5p(record: &bam::Record) -> i64 {
    if record.flags().is_reverse_complemented() {
        unclipped_end(record)
    } else {
        unclipped_start(record)
    }
}

fn tid(record: &bam::Record) -> i64 {
    record
        .reference_sequence_id()
        .and_then(std::result::Result::ok)
        .map_or(-1, |t| t as i64)
}

fn mate_tid(record: &bam::Record) -> i64 {
    record
        .mate_reference_sequence_id()
        .and_then(std::result::Result::ok)
        .map_or(-1, |t| t as i64)
}

fn base_qual_sum(record: &bam::Record) -> i64 {
    // samtools markdup sums only base qualities >= 15 (its default threshold).
    record
        .quality_scores()
        .as_ref()
        .iter()
        .filter(|&&q| q >= 15)
        .map(|&q| i64::from(q))
        .sum()
}

fn read_name(record: &bam::Record) -> Vec<u8> {
    record.name().map_or_else(Vec::new, |n| {
        let bytes: &[u8] = n.as_ref();
        bytes.to_vec()
    })
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

fn single_key_for(record: &bam::Record) -> SingleKey {
    SingleKey {
        tid: tid(record),
        coord: unclipped_5p(record),
        rev: record.flags().is_reverse_complemented(),
    }
}

/// Build the pair key from this read's perspective, given the mate's unclipped_5' coord.
/// Mirrors samtools template-mode make_pair_key (no barcode/RG).
fn pair_key_for(record: &bam::Record, mate_5p: i64) -> PairKey {
    let f = record.flags();
    let this_ref = tid(record);
    let other_ref = mate_tid(record);
    let this_rev = f.is_reverse_complemented();
    let mate_rev = f.is_mate_reverse_complemented();

    // Compute unclipped start and end for this read.
    let this_start = unclipped_start(record);
    let this_end = unclipped_end(record);
    // Mate's unclipped 5' is mate_5p. We don't have the mate's full span, so
    // we use mate_5p as both the "other_coord" and accept that this approximates
    // samtools' template mode which uses both other_start and other_end.
    // For the standard FR orientation (most common), samtools uses:
    //   forward read: this_coord = unclipped_start, other_coord = mate_unclipped_end
    //   reverse read: this_coord = unclipped_end,   other_coord = mate_unclipped_start
    // Since we have the mate record available (via our name-index), we derive
    // other_start = mate_5p for forward mate, other_end = mate_5p for reverse mate,
    // and vice versa for the other coordinate.  This matches unclipped_5p() semantics.

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

    // Orientation and coordinate selection (mirroring samtools template mode):
    //   leftmost + (!this_rev && mate_rev) → O_FR, this_coord=this_start, other_coord=mate_5p(=end)
    //   leftmost + (this_rev && !mate_rev) → O_RF, this_coord=this_end,   other_coord=mate_5p(=start)
    //   !leftmost + (!this_rev && mate_rev) → O_RF, this_coord=this_start, other_coord=mate_5p
    //   !leftmost + (this_rev && !mate_rev) → O_FR, this_coord=this_end,   other_coord=mate_5p
    // For same-orientation FF/RR the orientation encodes read1/read2 role; we simplify
    // to just (this_rev, mate_rev) which is sufficient to distinguish the four classes.
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
#[derive(Debug, Clone)]
struct SingleEntry {
    idx: usize,
    is_paired: bool,
}

pub fn markdup(
    input: &Path,
    output: &mut dyn Write,
    opts: &MarkdupOpts,
    workers: NonZero<usize>,
) -> Result<MarkdupStats> {
    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    let mut records: Vec<bam::Record> = Vec::new();
    for result in reader.records() {
        records.push(result.map_err(RsomicsError::Io)?);
    }

    // --- Pass 1: resolve mate unclipped-5' coordinates for paired reads. ---
    //
    // samtools markdup requires the MC (mate CIGAR) tag from fixmate to compute
    // the mate's unclipped position.  We instead do a name-grouping pass over
    // the already-loaded records, giving us both ends without requiring fixmate.
    let mut mate_5p_map: HashMap<usize, i64> = HashMap::new();
    {
        let mut name_to_idxs: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
        for (i, r) in records.iter().enumerate() {
            let f = r.flags();
            if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
                continue;
            }
            if f.is_segmented() && !f.is_mate_unmapped() {
                name_to_idxs.entry(read_name(r)).or_default().push(i);
            }
        }
        for idxs in name_to_idxs.values() {
            if idxs.len() != 2 {
                // Chimeric / multi-segment: treat each end as single-end.
                continue;
            }
            let (i0, i1) = (idxs[0], idxs[1]);
            mate_5p_map.insert(i0, unclipped_5p(&records[i1]));
            mate_5p_map.insert(i1, unclipped_5p(&records[i0]));
        }
    }

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
        if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
            continue;
        }

        let sk = single_key_for(record);
        let score = base_qual_sum(record);
        let is_paired = f.is_segmented() && !f.is_mate_unmapped() && mate_5p_map.contains_key(&i);

        // --- single_hash: every read registers its own-end key ---
        match single_hash.get(&sk).cloned() {
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
                let old_score = base_qual_sum(&records[existing.idx]);
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
        if is_paired {
            let mate_5p = mate_5p_map[&i];
            let pk = pair_key_for(record, mate_5p);

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

    let mut writer = bam::io::Writer::new(output);
    writer.write_header(&header).map_err(RsomicsError::Io)?;

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
            let mut buf =
                RecordBuf::try_from_alignment_record(&header, record).map_err(RsomicsError::Io)?;
            *buf.flags_mut() |= sam::alignment::record::Flags::DUPLICATE;
            writer
                .write_alignment_record(&header, &buf)
                .map_err(RsomicsError::Io)?;
        } else {
            writer
                .write_record(&header, record)
                .map_err(RsomicsError::Io)?;
        }
    }

    Ok(stats)
}
