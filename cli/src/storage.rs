use anyhow::Context;
use anyhow::Result;
use serde_json::Value;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command as ProcessCommand;
use std::process::Stdio;
use std::time::SystemTime;
use std::time::UNIX_EPOCH;

pub(super) fn now_unix_secs() -> Result<u64> {
    Ok(SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system time is before the Unix epoch")?
        .as_secs())
}

pub(super) fn read_json_file(path: &Path, label: &str) -> Result<Value> {
    serde_json::from_str(
        &fs::read_to_string(path)
            .with_context(|| format!("failed to read {label} {}", path.display()))?,
    )
    .with_context(|| format!("invalid JSON in {label} {}", path.display()))
}

/// Creates `path` exclusively (O_EXCL) and writes `bytes`. Refuses to follow a
/// pre-planted file or symlink at a predictable temporary name.
pub(super) fn write_new_file(path: &Path, bytes: &[u8]) -> Result<()> {
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    file.write_all(bytes)
        .with_context(|| format!("failed to write {}", path.display()))
}

/// Unique sibling name for an atomic write: pid, wall clock and a per-process
/// counter, so concurrent writers in one process never collide.
fn temporary_sibling(path: &Path, fallback: &str, suffix: &str) -> Result<PathBuf> {
    static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or(fallback);
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system time is before the Unix epoch")?
        .as_nanos();
    let sequence = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    Ok(path.with_file_name(format!(
        ".{name}.tmp-{}-{nanos}-{sequence}{suffix}",
        std::process::id()
    )))
}

pub(super) fn write_json_atomic(path: &Path, value: &Value) -> Result<()> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let temporary = temporary_sibling(path, "state.json", "")?;
    write_new_file(&temporary, &serde_json::to_vec_pretty(value)?)?;
    sync_file(&temporary)?;
    if let Err(error) = fs::rename(&temporary, path) {
        let _ = fs::remove_file(&temporary);
        return Err(error)
            .with_context(|| format!("failed to atomically replace {}", path.display()));
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        sync_directory(parent)?;
    }
    Ok(())
}

pub(super) fn write_bytes_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let temporary = temporary_sibling(path, "export", "")?;
    write_new_file(&temporary, bytes)?;
    sync_file(&temporary)?;
    if let Err(error) = fs::rename(&temporary, path) {
        let _ = fs::remove_file(&temporary);
        return Err(error)
            .with_context(|| format!("failed to atomically replace {}", path.display()));
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        sync_directory(parent)?;
    }
    Ok(())
}

pub(super) fn path_with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut value = path.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

pub(super) fn write_bgzip_vcf_with_tabix(path: &Path, vcf: &str) -> Result<PathBuf> {
    if path.extension().and_then(|value| value.to_str()) != Some("gz") {
        anyhow::bail!("--bgzip-vcf output must end in .gz: {}", path.display());
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    let temporary = temporary_sibling(path, "variants.vcf.gz", ".gz")?;
    let temporary_index = path_with_suffix(&temporary, ".tbi");
    let final_index = path_with_suffix(path, ".tbi");

    let mut child = ProcessCommand::new("bgzip")
        .arg("-c")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .context("failed to start bgzip; install htslib/bgzip")?;
    // bgzip streams compressed blocks while it reads, so feed stdin from a
    // separate thread: writing everything before draining stdout deadlocks
    // once the output exceeds the pipe buffer (about 64 KiB).
    let mut stdin = child.stdin.take().context("bgzip stdin is unavailable")?;
    let compressed = std::thread::scope(|scope| {
        let writer = scope.spawn(move || stdin.write_all(vcf.as_bytes()));
        let output = child.wait_with_output();
        let written = writer
            .join()
            .map_err(|_| anyhow::anyhow!("bgzip input writer panicked"))?;
        written.context("failed to stream VCF into bgzip")?;
        output.context("failed to wait for bgzip")
    })?;
    if !compressed.status.success() {
        anyhow::bail!(
            "bgzip failed: {}",
            String::from_utf8_lossy(&compressed.stderr).trim()
        );
    }
    write_new_file(&temporary, &compressed.stdout)?;
    sync_file(&temporary)?;

    let tabix = ProcessCommand::new("tabix")
        .args(["-f", "-p", "vcf"])
        .arg(&temporary)
        .output()
        .context("failed to start tabix; install htslib/tabix")?;
    if !tabix.status.success() {
        let _ = fs::remove_file(&temporary);
        let _ = fs::remove_file(&temporary_index);
        anyhow::bail!(
            "tabix failed: {}",
            String::from_utf8_lossy(&tabix.stderr).trim()
        );
    }
    sync_file(&temporary_index)?;
    if let Err(error) = fs::rename(&temporary, path) {
        let _ = fs::remove_file(&temporary);
        let _ = fs::remove_file(&temporary_index);
        return Err(error).with_context(|| format!("failed to publish {}", path.display()));
    }
    if let Err(error) = fs::rename(&temporary_index, &final_index) {
        let _ = fs::remove_file(path);
        let _ = fs::remove_file(&temporary_index);
        return Err(error).with_context(|| format!("failed to publish {}", final_index.display()));
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        sync_directory(parent)?;
    }
    Ok(final_index)
}

pub(super) fn sha256_bytes(bytes: &[u8]) -> Result<String> {
    let mut child = ProcessCommand::new("sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .context("failed to start sha256sum")?;
    child
        .stdin
        .as_mut()
        .context("sha256sum stdin is unavailable")?
        .write_all(bytes)
        .context("failed to stream bytes into sha256sum")?;
    let output = child
        .wait_with_output()
        .context("failed to wait for sha256sum")?;
    if !output.status.success() {
        anyhow::bail!("sha256sum failed while hashing analysis result");
    }
    String::from_utf8(output.stdout)
        .context("sha256sum returned non-UTF-8 output")?
        .split_whitespace()
        .next()
        .map(str::to_owned)
        .filter(|digest| digest.len() == 64 && digest.bytes().all(|byte| byte.is_ascii_hexdigit()))
        .context("sha256sum returned an invalid digest")
}

pub(super) fn sync_file(path: &Path) -> Result<()> {
    fs::File::open(path)
        .with_context(|| format!("failed to open {} for synchronization", path.display()))?
        .sync_all()
        .with_context(|| format!("failed to synchronize {}", path.display()))
}

pub(super) fn sync_directory(path: &Path) -> Result<()> {
    fs::File::open(path)
        .with_context(|| {
            format!(
                "failed to open directory {} for synchronization",
                path.display()
            )
        })?
        .sync_all()
        .with_context(|| format!("failed to synchronize directory {}", path.display()))
}

pub(super) fn normalize_text(value: &str) -> String {
    value.replace("\r\n", "\n").trim_end().to_string()
}

pub(super) fn sha256_file(path: &Path) -> Result<String> {
    let output = ProcessCommand::new("sha256sum")
        .arg(path)
        .output()
        .with_context(|| format!("failed to start sha256sum for {}", path.display()))?;
    if !output.status.success() {
        anyhow::bail!("sha256sum failed for {}", path.display());
    }
    let text = String::from_utf8(output.stdout).context("sha256sum returned non-UTF-8 output")?;
    text.split_whitespace()
        .next()
        .filter(|value| value.len() == 64 && value.chars().all(|c| c.is_ascii_hexdigit()))
        .map(str::to_string)
        .context("sha256sum returned an invalid digest")
}
