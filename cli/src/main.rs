mod analyze;
mod architecture_contract;
mod args;
mod clinical;
mod comparator;
mod contract;
mod doctor;
mod exports;
mod fixtures;
mod haplogroup_validation;
mod projection;
mod report;
mod storage;
#[cfg(test)]
mod tests;

use crate::analyze::analyze;
use crate::args::{Cli, Command};
use crate::clinical::{manage_clinical_snapshot, update_clinical};
use crate::comparator::compare_mtdna_server2;
use crate::doctor::doctor;
use crate::fixtures::{
    validate_clinical_fixture, validate_clinical_manifest, validate_error_manifest,
    validate_evidence_fixture, validate_evidence_graph_fixture, validate_fixture,
    validate_sv_fixture,
};
use crate::haplogroup_validation::validate_haplogroup_manifest;
use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Analyze(args) => analyze(*args),
        Command::CompareMtDnaServer2(args) => compare_mtdna_server2(args),
        Command::Doctor => doctor(),
        Command::UpdateClinical(args) => update_clinical(args).map(|path| {
            println!("clinical annotation cache: {}", path.display());
            println!("set MITO_CLINICAL_ANNOTATIONS={} to use it", path.display());
        }),
        Command::ClinicalSnapshot(args) => manage_clinical_snapshot(args),
        Command::ValidateFixture(args) => validate_fixture(args),
        Command::ValidateClinicalFixture(args) => validate_clinical_fixture(args),
        Command::ValidateClinicalManifest(args) => validate_clinical_manifest(args),
        Command::ValidateEvidenceFixture(args) => validate_evidence_fixture(args),
        Command::ValidateEvidenceGraphFixture(args) => validate_evidence_graph_fixture(args),
        Command::ValidateErrorManifest(args) => validate_error_manifest(args),
        Command::ValidateHaplogroupManifest(args) => validate_haplogroup_manifest(args),
        Command::ValidateSvFixture(args) => validate_sv_fixture(args),
    }
}
