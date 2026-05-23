use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam-markdup"))
}

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn samtools_available() -> bool {
    Command::new("samtools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn run_ok(cmd: &mut Command) {
    assert!(cmd.status().unwrap().success(), "command failed: {cmd:?}");
}

/// Sorted (QNAME, FLAG) pairs of reads flagged as duplicates (0x400).
fn dup_records(bam: &Path) -> Vec<(String, u16)> {
    let out = Command::new("samtools")
        .args(["view", "-f", "1024"])
        .arg(bam)
        .output()
        .unwrap();
    assert!(out.status.success());
    let mut records: Vec<(String, u16)> = String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter_map(|l| {
            let mut cols = l.split('\t');
            let name = cols.next()?.to_owned();
            let flag: u16 = cols.next()?.parse().ok()?;
            Some((name, flag))
        })
        .collect();
    records.sort();
    records
}

/// Run samtools fixmate+markdup pipeline and our tool on the same SAM, then
/// compare the sets of reads flagged as duplicates.
fn run_compat(sam: &Path, dir: &Path, tag: &str) {
    let all = dir.join(format!("{tag}_all.bam"));
    {
        let f = std::fs::File::create(&all).unwrap();
        assert!(
            Command::new("samtools")
                .args(["view", "-b"])
                .arg(sam)
                .stdout(f)
                .status()
                .unwrap()
                .success()
        );
    }
    // samtools pipeline: name-sort | fixmate -m | coord-sort | markdup
    let ns = dir.join(format!("{tag}_ns.bam"));
    run_ok(
        Command::new("samtools")
            .args(["sort", "-n", "-o"])
            .arg(&ns)
            .arg(&all),
    );
    let fm = dir.join(format!("{tag}_fm.bam"));
    run_ok(
        Command::new("samtools")
            .args(["fixmate", "-m"])
            .arg(&ns)
            .arg(&fm),
    );
    let cs = dir.join(format!("{tag}_cs.bam"));
    run_ok(
        Command::new("samtools")
            .args(["sort", "-o"])
            .arg(&cs)
            .arg(&fm),
    );
    let smd = dir.join(format!("{tag}_smd.bam"));
    run_ok(Command::new("samtools").arg("markdup").arg(&cs).arg(&smd));

    // our tool operates on the same coordinate-sorted BAM (no fixmate needed)
    let omd = dir.join(format!("{tag}_omd.bam"));
    run_ok(ours().arg(&cs).arg("-o").arg(&omd));

    assert_eq!(
        dup_records(&omd),
        dup_records(&smd),
        "[{tag}] duplicate-flagged read set must match samtools markdup"
    );
}

// Pure paired-end: rA is original, rD has same unclipped 5' due to soft-clip.
// Exercises unclipped-coordinate logic on all-PE data.
#[test]
fn markdup_pure_pe_matches_samtools() {
    if !samtools_available() {
        eprintln!("skipping: samtools not found");
        return;
    }
    let dir = std::env::temp_dir().join("rsomics-bam-markdup-compat-pe");
    let _ = std::fs::create_dir_all(&dir);
    run_compat(&golden_dir().join("small_pe.sam"), &dir, "pe");
    let _ = std::fs::remove_dir_all(&dir);
}

// Mixed SE+PE: rA is an original PE pair; rB is a duplicate PE pair; rC is a
// SE read at the same position as rA (PE beats SE); rE/rF are two SE reads at
// the same position (higher-quality rE wins).  Expected dups: rB/1, rB/2, rC, rF.
#[test]
fn markdup_mixed_se_pe_matches_samtools() {
    if !samtools_available() {
        eprintln!("skipping: samtools not found");
        return;
    }
    let dir = std::env::temp_dir().join("rsomics-bam-markdup-compat-mixed");
    let _ = std::fs::create_dir_all(&dir);
    run_compat(&golden_dir().join("mixed_se_pe.sam"), &dir, "mixed");
    let _ = std::fs::remove_dir_all(&dir);
}
