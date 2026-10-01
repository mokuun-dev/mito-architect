use anyhow::Result;
use mito_ffi::MitoEngine;
use std::process::Command as ProcessCommand;

pub(super) fn doctor() -> Result<()> {
    let capabilities = MitoEngine::capabilities();
    println!(
        "mito-engine {} (result schema {}, error schema {})",
        capabilities.engine_version, capabilities.schema_version, capabilities.error_schema_version
    );
    println!(
        "BAM/CRAM reader: {}",
        if capabilities.htslib {
            "enabled"
        } else {
            "disabled"
        }
    );

    let mut missing = Vec::new();
    for (program, arguments, required) in [
        ("samtools", &["--version"][..], true),
        ("minimap2", &["--version"][..], true),
        ("bcftools", &["--version"][..], false),
        ("bgzip", &["--version"][..], false),
        ("tabix", &["--version"][..], false),
        ("fasterq-dump", &["--version"][..], false),
    ] {
        match ProcessCommand::new(program).args(arguments).output() {
            Ok(output) if output.status.success() => {
                let version = String::from_utf8_lossy(&output.stdout)
                    .lines()
                    .find(|line| !line.trim().is_empty())
                    .unwrap_or("version unavailable")
                    .trim()
                    .to_owned();
                println!("{program}: {version}");
            }
            _ => {
                println!("{program}: missing");
                if required {
                    missing.push(program);
                }
            }
        }
    }

    if !capabilities.htslib {
        anyhow::bail!(
            "native core was built without htslib; clean and rebuild after installing htslib"
        );
    }
    if !missing.is_empty() {
        anyhow::bail!("missing required tools: {}", missing.join(", "));
    }
    println!("required native analysis capabilities are available");
    Ok(())
}
