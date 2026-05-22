use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam-markdup"))
}

fn golden() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small_pe.sam")
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

/// Sorted QNAMEs of reads flagged as duplicates (0x400).
fn dup_names(bam: &Path) -> Vec<String> {
    let out = Command::new("samtools")
        .args(["view", "-f", "1024"])
        .arg(bam)
        .output()
        .unwrap();
    assert!(out.status.success());
    let mut names: Vec<String> = String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter_map(|l| l.split('\t').next().map(str::to_owned))
        .collect();
    names.sort();
    names
}

// ours must flag the same reads as `samtools markdup`. The golden has a plain
// duplicate pair (rB) and a soft-clipped duplicate pair (rD, where clipped pos
// differs but unclipped 5' matches) — exercising the unclipped-coordinate logic.
#[test]
fn markdup_matches_samtools() {
    if !samtools_available() {
        eprintln!("skipping: samtools not found");
        return;
    }
    let dir = std::env::temp_dir().join("rsomics-bam-markdup-compat");
    let _ = std::fs::create_dir_all(&dir);

    let all = dir.join("all.bam");
    {
        let f = std::fs::File::create(&all).unwrap();
        assert!(
            Command::new("samtools")
                .args(["view", "-b"])
                .arg(golden())
                .stdout(f)
                .status()
                .unwrap()
                .success()
        );
    }
    // samtools markdup pipeline: name-sort | fixmate -m | coord-sort | markdup
    let ns = dir.join("ns.bam");
    run_ok(
        Command::new("samtools")
            .args(["sort", "-n", "-o"])
            .arg(&ns)
            .arg(&all),
    );
    let fm = dir.join("fm.bam");
    run_ok(
        Command::new("samtools")
            .args(["fixmate", "-m"])
            .arg(&ns)
            .arg(&fm),
    );
    let cs = dir.join("cs.bam");
    run_ok(
        Command::new("samtools")
            .args(["sort", "-o"])
            .arg(&cs)
            .arg(&fm),
    );
    let smd = dir.join("smd.bam");
    run_ok(Command::new("samtools").arg("markdup").arg(&cs).arg(&smd));

    // ours operates on the same coordinate-sorted BAM (no fixmate needed)
    let omd = dir.join("omd.bam");
    run_ok(ours().arg(&cs).arg("-o").arg(&omd));

    assert_eq!(
        dup_names(&omd),
        dup_names(&smd),
        "duplicate-flagged read set must match samtools markdup"
    );
}
