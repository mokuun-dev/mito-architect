use clap::Args;
use clap::Parser;
use clap::Subcommand;
use clap::ValueEnum;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "mito-cli")]
#[command(about = "Long-read mtDNA analysis and offline report generator")]
pub(super) struct Cli {
    #[command(subcommand)]
    pub(super) command: Command,
}

#[derive(Debug, Subcommand)]
pub(super) enum Command {
    Analyze(Box<AnalyzeArgs>),
    CompareMtDnaServer2(CompareMtDnaServer2Args),
    Doctor,
    UpdateClinical(UpdateClinicalArgs),
    ClinicalSnapshot(ClinicalSnapshotArgs),
    ValidateFixture(ValidateFixtureArgs),
    ValidateClinicalFixture(ValidateClinicalFixtureArgs),
    ValidateClinicalManifest(ValidateClinicalManifestArgs),
    ValidateEvidenceFixture(ValidateEvidenceFixtureArgs),
    ValidateEvidenceGraphFixture(ValidateEvidenceGraphFixtureArgs),
    ValidateErrorManifest(ValidateErrorManifestArgs),
    ValidateHaplogroupManifest(ValidateHaplogroupManifestArgs),
    ValidateSvFixture(ValidateSvFixtureArgs),
}

#[derive(Clone, Copy, Debug, ValueEnum)]
pub(super) enum ComparatorFilterPolicy {
    PassOnly,
    All,
}

#[derive(Debug, Args)]
pub(super) struct CompareMtDnaServer2Args {
    /// Authoritative Mito-Architect schema 0.6 JSON result.
    #[arg(long = "result")]
    pub(super) result: PathBuf,

    /// mtDNA-Server 2 variants.annotated.txt file.
    #[arg(long = "comparator")]
    pub(super) comparator: PathBuf,

    /// Exact mtDNA-Server 2 sample ID. Required for multi-sample files.
    #[arg(long = "sample")]
    pub(super) sample: Option<String>,

    /// Pinned upstream release whose variants.annotated.txt format is parsed.
    #[arg(long = "comparator-version", default_value = "2.1.16")]
    pub(super) comparator_version: String,

    /// Include only PASS comparator records or every filter state.
    #[arg(long = "filter-policy", value_enum, default_value_t = ComparatorFilterPolicy::PassOnly)]
    pub(super) filter_policy: ComparatorFilterPolicy,

    /// Deterministic machine-readable differential report.
    #[arg(short = 'o', long = "output")]
    pub(super) output: PathBuf,

    /// Optional development gate; comparator concordance is not truth accuracy.
    #[arg(long = "min-call-concordance")]
    pub(super) min_call_concordance: Option<f64>,

    /// Optional development gate over variants called by both tools.
    #[arg(long = "max-mean-hf-delta")]
    pub(super) max_mean_hf_delta: Option<f64>,
}

#[derive(Debug, Args)]
pub(super) struct AnalyzeArgs {
    #[arg(short = 'i', long = "input")]
    pub(super) input: PathBuf,

    #[arg(short = 'o', long = "output", default_value = "output.html")]
    pub(super) output: PathBuf,

    #[arg(short = 'r', long = "reference")]
    pub(super) reference: Option<PathBuf>,

    #[arg(long = "json")]
    pub(super) json: bool,

    #[arg(long = "vcf")]
    pub(super) vcf: Option<PathBuf>,

    /// Write the unified schema 0.6 SNV/small-indel projection as TSV.
    #[arg(long = "tsv")]
    pub(super) tsv: Option<PathBuf>,

    /// Write bgzip-compressed VCF and a tabix index (requires bgzip and tabix).
    #[arg(long = "bgzip-vcf")]
    pub(super) bgzip_vcf: Option<PathBuf>,

    /// Write a deterministic per-analysis provenance manifest.
    #[arg(long = "provenance-manifest")]
    pub(super) provenance_manifest: Option<PathBuf>,

    /// Optional directory for deterministic schema 0.6 observation-page sidecars.
    #[arg(long = "evidence-pages-dir")]
    pub(super) evidence_pages_dir: Option<PathBuf>,

    #[arg(long = "filter-numt", default_value_t = true)]
    pub(super) filter_numt: bool,

    #[arg(long = "threads", default_value_t = 1)]
    pub(super) threads: usize,

    #[arg(long = "min-mapq", default_value_t = 20)]
    pub(super) min_mapping_quality: u8,

    #[arg(long = "min-base-quality", default_value_t = 10)]
    pub(super) min_base_quality: u8,

    #[arg(long = "excluded-snp-flags", default_value_t = 3840)]
    pub(super) excluded_snp_flags: u16,

    #[arg(long = "numt-threshold", default_value_t = 0.30)]
    pub(super) numt_threshold: f64,

    #[arg(long = "allow-development-tags", default_value_t = false)]
    pub(super) allow_development_tags: bool,

    /// Emit the opt-in schema 0.6 fragment/molecule/event evidence graph.
    #[arg(long = "evidence-graph", default_value_t = false)]
    pub(super) emit_evidence_graph: bool,

    /// Hard cap for in-memory schema 0.6 sparse observations.
    #[arg(long = "max-evidence-observations", default_value_t = 5_000_000)]
    pub(super) max_evidence_observations: usize,

    /// Hard cap for callable-aware event-pair phase projections.
    #[arg(long = "max-phase-links", default_value_t = 1_000_000)]
    pub(super) max_phase_links: usize,

    /// Hard cap for pair evaluations while building the schema 0.6 phase graph.
    #[arg(long = "max-phase-work", default_value_t = 20_000_000)]
    pub(super) max_phase_work: usize,

    /// Hard cap for molecule references retained by phase-link evidence.
    #[arg(long = "max-phase-molecule-references", default_value_t = 5_000_000)]
    pub(super) max_phase_molecule_references: usize,

    /// Hard cap for the serialized analysis JSON payload.
    #[arg(long = "max-result-bytes", default_value_t = 128 * 1024 * 1024)]
    pub(super) max_result_bytes: usize,

    /// Maximum observations per schema 0.6 columnar page.
    #[arg(long = "evidence-page-size", default_value_t = 4096)]
    pub(super) evidence_page_size: usize,

    /// Minimum independent molecules required to seed and retain a candidate architecture.
    #[arg(long = "min-architecture-molecules", default_value_t = 2)]
    pub(super) min_architecture_molecules: usize,

    /// Hard cap for callable-aware candidate architectures.
    #[arg(long = "max-candidate-architectures", default_value_t = 256)]
    pub(super) max_candidate_architectures: usize,

    /// Maximum callable-aware disagreement fraction for assignment.
    #[arg(long = "architecture-max-distance", default_value_t = 0.20)]
    pub(super) architecture_max_distance: f64,

    /// Best-versus-second score margin reported as ambiguous.
    #[arg(long = "architecture-ambiguity-margin", default_value_t = 0.05)]
    pub(super) architecture_ambiguity_margin: f64,

    /// Minimum weighted candidate profile overlap required for assignment.
    #[arg(long = "architecture-min-overlap", default_value_t = 0.50)]
    pub(super) architecture_min_overlap_fraction: f64,

    /// Alternate fraction required for a defining event.
    #[arg(long = "architecture-consensus", default_value_t = 0.80)]
    pub(super) architecture_consensus_fraction: f64,

    /// Alternate fraction required for an optional event.
    #[arg(long = "architecture-optional", default_value_t = 0.20)]
    pub(super) architecture_optional_fraction: f64,

    /// Deterministic signature-stability resamples; zero disables the estimate.
    #[arg(long = "architecture-stability-replicates", default_value_t = 32)]
    pub(super) architecture_stability_replicates: usize,

    /// Recorded deterministic stability seed.
    #[arg(
        long = "architecture-seed",
        default_value_t = 0x4d49_544f_4152_4348_u64
    )]
    pub(super) architecture_seed: u64,

    /// Explicit SAM tag used as the physical-molecule identifier (for example MI).
    #[arg(long = "molecule-id-tag", default_value = "")]
    pub(super) molecule_id_tag: String,

    /// Optional SAM tag carrying UMI metadata (for example RX).
    #[arg(long = "umi-tag", default_value = "")]
    pub(super) umi_tag: String,

    /// Optional SAM tag carrying duplex metadata.
    #[arg(long = "duplex-tag", default_value = "")]
    pub(super) duplex_tag: String,

    #[arg(long = "update-clinical")]
    pub(super) update_clinical: bool,
}

#[derive(Debug, Args)]
pub(super) struct UpdateClinicalArgs {
    #[arg(long = "source")]
    pub(super) source: Option<PathBuf>,

    #[arg(long = "output")]
    pub(super) output: Option<PathBuf>,

    #[arg(long = "clinvar-live")]
    pub(super) clinvar_live: bool,

    #[arg(long = "clinvar-gz")]
    pub(super) clinvar_gz: Option<PathBuf>,
}

#[derive(Clone, Debug, ValueEnum)]
pub(super) enum ClinicalSnapshotAction {
    Stage,
    Activate,
    Rollback,
    Verify,
    Status,
}

#[derive(Debug, Args)]
pub(super) struct ClinicalSnapshotArgs {
    #[arg(value_enum)]
    pub(super) action: ClinicalSnapshotAction,

    #[arg(long = "store")]
    pub(super) store: Option<PathBuf>,

    #[arg(long = "source")]
    pub(super) source: Option<PathBuf>,

    #[arg(long = "snapshot-id")]
    pub(super) snapshot_id: Option<String>,

    #[arg(long = "license-id")]
    pub(super) license_id: Option<String>,

    #[arg(long = "source-policy")]
    pub(super) source_policy: Option<String>,

    #[arg(long = "activate")]
    pub(super) activate: bool,

    #[arg(long = "max-age-days")]
    pub(super) max_age_days: Option<u64>,
}

#[derive(Debug, Args)]
pub(super) struct ValidateFixtureArgs {
    #[arg(long = "input")]
    pub(super) input: PathBuf,

    #[arg(long = "expected-vcf")]
    pub(super) expected_vcf: PathBuf,

    #[arg(long = "expected-passed")]
    pub(super) expected_passed: u64,

    #[arg(long = "expected-numt")]
    pub(super) expected_numt: u64,

    #[arg(long = "expected-snp", required = true, action = clap::ArgAction::Append)]
    pub(super) expected_snp: Vec<String>,

    #[arg(long = "expected-sv", action = clap::ArgAction::Append)]
    pub(super) expected_sv: Vec<String>,

    #[arg(long = "expected-mapq", action = clap::ArgAction::Append)]
    pub(super) expected_mapq: Vec<u8>,

    #[arg(long = "expected-aux", action = clap::ArgAction::Append)]
    pub(super) expected_aux: Vec<String>,
}

#[derive(Debug, Args)]
pub(super) struct ValidateSvFixtureArgs {
    #[arg(long = "input")]
    pub(super) input: PathBuf,

    #[arg(long = "expected-json")]
    pub(super) expected_json: PathBuf,
}

#[derive(Debug, Args)]
pub(super) struct ValidateEvidenceFixtureArgs {
    #[arg(long = "input")]
    pub(super) input: PathBuf,

    #[arg(long = "expected-json")]
    pub(super) expected_json: PathBuf,

    /// Write the deterministic projection before comparison. Intended for
    /// reviewed fixture regeneration after an explicit contract change.
    #[arg(long = "write-projection")]
    pub(super) write_projection: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub(super) struct ValidateEvidenceGraphFixtureArgs {
    #[arg(long = "input")]
    pub(super) input: PathBuf,

    #[arg(long = "expected-json")]
    pub(super) expected_json: PathBuf,

    #[arg(long = "evidence-page-size", default_value_t = 3)]
    pub(super) evidence_page_size: usize,

    #[arg(long = "molecule-id-tag", default_value = "")]
    pub(super) molecule_id_tag: String,

    #[arg(long = "umi-tag", default_value = "")]
    pub(super) umi_tag: String,

    #[arg(long = "duplex-tag", default_value = "")]
    pub(super) duplex_tag: String,

    /// Write the deterministic projection before comparison. Intended for
    /// reviewed fixture regeneration, not normal analysis output.
    #[arg(long = "write-projection")]
    pub(super) write_projection: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub(super) struct ValidateClinicalFixtureArgs {
    #[arg(long = "input")]
    pub(super) input: PathBuf,

    #[arg(long = "annotations")]
    pub(super) annotations: PathBuf,

    #[arg(long = "expected-json")]
    pub(super) expected_json: PathBuf,
}

#[derive(Debug, Args)]
pub(super) struct ValidateClinicalManifestArgs {
    #[arg(long = "manifest")]
    pub(super) manifest: PathBuf,
}

#[derive(Debug, Args)]
pub(super) struct ValidateErrorManifestArgs {
    #[arg(long = "manifest")]
    pub(super) manifest: PathBuf,
}

#[derive(Debug, Args)]
pub(super) struct ValidateHaplogroupManifestArgs {
    #[arg(long = "manifest")]
    pub(super) manifest: PathBuf,
}
