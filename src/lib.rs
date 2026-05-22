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

/// Unclipped 5' coordinate: forward = pos - leading clip; reverse = pos + ref_span + trailing clip.
fn unclipped_5p(record: &bam::Record) -> i64 {
    let pos = record
        .alignment_start()
        .and_then(std::result::Result::ok)
        .map_or(0, |p| p.get() as i64);
    let (span, lead, trail) = span_and_clips(record);
    if record.flags().is_reverse_complemented() {
        pos + span + trail
    } else {
        pos - lead
    }
}

fn tid(record: &bam::Record) -> i64 {
    record
        .reference_sequence_id()
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

fn name_key(record: &bam::Record) -> Vec<u8> {
    record.name().map_or_else(Vec::new, |n| {
        let bytes: &[u8] = n.as_ref();
        bytes.to_vec()
    })
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

    // group primary mapped reads by read name
    let mut groups: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
    for (i, r) in records.iter().enumerate() {
        let f = r.flags();
        if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
            continue;
        }
        groups.entry(name_key(r)).or_default().push(i);
    }

    // signature -> (best score, best group's member indices)
    let mut best: HashMap<String, (i64, Vec<usize>)> = HashMap::new();
    let mut is_dup = vec![false; records.len()];

    // Process groups in coordinate order (records are coordinate-sorted, so the
    // smallest member index = earliest coordinate). On equal score this keeps the
    // first-by-coordinate representative, matching samtools markdup's tie-break.
    let mut group_list: Vec<&Vec<usize>> = groups.values().collect();
    group_list.sort_by_key(|m| m.iter().min().copied().unwrap_or(usize::MAX));

    for members in group_list {
        // signature: for a pair, both ends' (tid, unclipped-5') sorted + orientation;
        // for a singleton, its (tid, unclipped-5', strand).
        let mut ends: Vec<(i64, i64, bool)> = members
            .iter()
            .map(|&i| {
                let r = &records[i];
                (tid(r), unclipped_5p(r), r.flags().is_reverse_complemented())
            })
            .collect();
        ends.sort_by_key(|e| (e.0, e.1));
        let sig = ends
            .iter()
            .map(|(t, u, rev)| format!("{t}:{u}:{}", u8::from(*rev)))
            .collect::<Vec<_>>()
            .join("|");
        let score: i64 = members.iter().map(|&i| base_qual_sum(&records[i])).sum();

        match best.get_mut(&sig) {
            Some(entry) => {
                if score > entry.0 {
                    for &i in &entry.1 {
                        is_dup[i] = true;
                    }
                    *entry = (score, members.clone());
                } else {
                    for &i in members {
                        is_dup[i] = true;
                    }
                }
            }
            None => {
                best.insert(sig, (score, members.clone()));
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
