//! Duplicate-detection keys, mirroring samtools `bam_markdup.c` (MIT) template
//! mode. The single key fingerprints one read by its unclipped 5' coordinate +
//! orientation; the pair key fingerprints a template by both ends' unclipped
//! coordinates, the leftmost flag and the pair orientation. Two reads of the
//! same template produce different pair keys (leftmost differs), so only reads
//! from distinct templates at identical positions collide — the duplicate
//! condition.

use rsomics_bamio::raw::RawRecord;

const FLAG_REVERSE: u16 = 0x10;
const FLAG_MATE_REVERSE: u16 = 0x20;
const FLAG_READ1: u16 = 0x40;

// CIGAR op codes (BAM packed encoding, low nibble): M=0 I=1 D=2 N=3 S=4 H=5 P=6 ==7 X=8.
const CIGAR_MATCH: u8 = 0;
const CIGAR_DELETION: u8 = 2;
const CIGAR_SKIP: u8 = 3;
const CIGAR_SOFT_CLIP: u8 = 4;
const CIGAR_HARD_CLIP: u8 = 5;
const CIGAR_SEQ_MATCH: u8 = 7;
const CIGAR_SEQ_MISMATCH: u8 = 8;

// Pair orientation codes (samtools O_FF/O_RR/O_FR/O_RF — the prime values are a
// hashing detail upstream; here only equality matters, so we keep the codes).
const O_FF: i8 = 2;
const O_RR: i8 = 3;
const O_FR: i8 = 5;
const O_RF: i8 = 7;

// Leftmost vs rightmost end (samtools R_LE/R_RI).
const R_LE: i8 = 11;
const R_RI: i8 = 13;

#[inline]
pub fn is_reverse(r: &RawRecord) -> bool {
    r.flags() & FLAG_REVERSE != 0
}

#[inline]
fn is_mate_reverse(r: &RawRecord) -> bool {
    r.flags() & FLAG_MATE_REVERSE != 0
}

#[inline]
fn is_read1(r: &RawRecord) -> bool {
    r.flags() & FLAG_READ1 != 0
}

/// Whether a CIGAR op consumes reference bases (M/D/N/=/X).
#[inline]
fn consumes_ref(op: u8) -> bool {
    matches!(
        op,
        CIGAR_MATCH | CIGAR_DELETION | CIGAR_SKIP | CIGAR_SEQ_MATCH | CIGAR_SEQ_MISMATCH
    )
}

/// Unclipped start (samtools `unclipped_start`): `pos - leading_clip + 1`,
/// 1-based. Both soft (S) and hard (H) clips count, scanning from the start
/// until the first non-clip op.
pub fn unclipped_start(r: &RawRecord) -> i64 {
    let mut clipped = 0i64;
    for (op, len) in r.cigar_ops() {
        if op == CIGAR_SOFT_CLIP || op == CIGAR_HARD_CLIP {
            clipped += i64::from(len);
        } else {
            break;
        }
    }
    i64::from(r.alignment_start()) - clipped + 1
}

/// Unclipped end (samtools `unclipped_end`): `bam_endpos(b) + trailing_clip`,
/// where `bam_endpos = pos + reference_span`. Trailing soft/hard clips count,
/// scanning from the end until the first non-clip op.
pub fn unclipped_end(r: &RawRecord) -> i64 {
    let mut span = 0i64;
    let mut trailing = 0i64;
    let ops: Vec<(u8, u32)> = r.cigar_ops().collect();
    for &(op, len) in &ops {
        if consumes_ref(op) {
            span += i64::from(len);
        }
    }
    for &(op, len) in ops.iter().rev() {
        if op == CIGAR_SOFT_CLIP || op == CIGAR_HARD_CLIP {
            trailing += i64::from(len);
        } else {
            break;
        }
    }
    i64::from(r.alignment_start()) + span + trailing
}

/// Mate unclipped start from the MC-tag CIGAR string and the mate's 0-based pos
/// (samtools `unclipped_other_start`): `mpos - leading_clip + 1`.
pub fn unclipped_other_start(mate_pos: i64, mc: &[u8]) -> i64 {
    let mut clipped = 0i64;
    let mut chars = mc.iter().peekable();
    while let Some(&&c) = chars.peek() {
        if c == b'*' || c == 0 {
            break;
        }
        let mut num: i64 = 0;
        if c.is_ascii_digit() {
            while let Some(&&d) = chars.peek() {
                if d.is_ascii_digit() {
                    num = num * 10 + i64::from(d - b'0');
                    chars.next();
                } else {
                    break;
                }
            }
        } else {
            num = 1;
        }
        match chars.next() {
            Some(b'S') | Some(b'H') => clipped += num,
            _ => break,
        }
    }
    mate_pos - clipped + 1
}

/// Mate unclipped end from the MC-tag CIGAR and the mate's 0-based pos
/// (samtools `unclipped_other_end`): `mpos + ref_span_with_trailing_clips`.
/// Reference-consuming ops add to the span; soft/hard clips add only once a
/// non-clip op has been seen (leading clips ignored).
pub fn unclipped_other_end(mate_pos: i64, mc: &[u8]) -> i64 {
    let mut refpos = 0i64;
    let mut skip = true; // ignore leading clips
    let mut chars = mc.iter().peekable();
    while let Some(&&c) = chars.peek() {
        if c == b'*' || c == 0 {
            break;
        }
        let mut num: i64 = 0;
        if c.is_ascii_digit() {
            while let Some(&&d) = chars.peek() {
                if d.is_ascii_digit() {
                    num = num * 10 + i64::from(d - b'0');
                    chars.next();
                } else {
                    break;
                }
            }
        } else {
            num = 1;
        }
        match chars.next() {
            Some(b'M') | Some(b'D') | Some(b'N') | Some(b'=') | Some(b'X') => {
                refpos += num;
                skip = false;
            }
            // Trailing clips count toward the unclipped end; leading clips (before
            // any reference-consuming op) are skipped, matching samtools.
            Some(b'S') | Some(b'H') if !skip => refpos += num,
            _ => {}
        }
    }
    mate_pos + refpos
}

/// Single-read key: `(this_ref, this_coord, orientation)` plus the per-read
/// `single` discriminant. `this_ref` is `tid + 1` (samtools offsets to keep 0
/// out of the hash; equality is preserved either way). Reverse reads key on the
/// unclipped end with O_RR; forward reads on the unclipped start with O_FF.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SingleKey {
    pub this_ref: i32,
    pub this_coord: i64,
    pub orientation: i8,
}

/// Pair key: both ends' unclipped coordinates, references, the leftmost flag and
/// the pair orientation. Built from one read's perspective via the MC tag, so
/// the two mates of a template yield distinct keys.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PairKey {
    pub this_ref: i32,
    pub this_coord: i64,
    pub other_ref: i32,
    pub other_coord: i64,
    pub leftmost: i8,
    pub orientation: i8,
}

/// Build the single key (samtools `make_single_key`). Returns the key and the
/// `this_coord` value samtools stores back as the read's window position.
pub fn make_single_key(r: &RawRecord) -> SingleKey {
    let this_ref = r.reference_sequence_id() + 1;
    let (this_coord, orientation) = if is_reverse(r) {
        (unclipped_end(r), O_RR)
    } else {
        (unclipped_start(r), O_FF)
    };
    SingleKey {
        this_ref,
        this_coord,
        orientation,
    }
}

/// Build the pair key in template mode (samtools `make_pair_key`,
/// `MD_MODE_TEMPLATE`). The MC tag supplies the mate CIGAR for the mate's
/// unclipped coordinates. Returns the key.
pub fn make_pair_key(r: &RawRecord, mc: &[u8]) -> PairKey {
    let this_ref = r.reference_sequence_id() + 1;
    let other_ref = r.mate_reference_sequence_id() + 1;

    let mut this_coord = unclipped_start(r);
    let this_end = unclipped_end(r);

    let mate_pos = i64::from(r.mate_alignment_start());
    let mut other_coord = unclipped_other_start(mate_pos, mc);
    let other_end = unclipped_other_end(mate_pos, mc);

    let this_rev = is_reverse(r);
    let mate_rev = is_mate_reverse(r);
    let read1 = is_read1(r);

    let leftmost: bool = if this_ref != other_ref {
        this_ref < other_ref
    } else if this_rev == mate_rev {
        if !this_rev {
            this_coord <= other_coord
        } else {
            this_end <= other_end
        }
    } else if this_rev {
        this_end <= other_coord
    } else {
        this_coord <= other_end
    };

    let orientation: i8 = if leftmost {
        if this_rev == mate_rev {
            other_coord = other_end;
            if !this_rev {
                if read1 { O_FF } else { O_RR }
            } else if read1 {
                O_RR
            } else {
                O_FF
            }
        } else if !this_rev {
            other_coord = other_end;
            O_FR
        } else {
            this_coord = this_end;
            O_RF
        }
    } else if this_rev == mate_rev {
        this_coord = this_end;
        if !this_rev {
            if read1 { O_RR } else { O_FF }
        } else if read1 {
            O_FF
        } else {
            O_RR
        }
    } else if !this_rev {
        other_coord = other_end;
        O_RF
    } else {
        this_coord = this_end;
        O_FR
    };

    let left_read = if leftmost { R_LE } else { R_RI };

    PairKey {
        this_ref,
        this_coord,
        other_ref,
        other_coord,
        leftmost: left_read,
        orientation,
    }
}
