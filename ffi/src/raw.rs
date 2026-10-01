use std::os::raw::{c_char, c_void};

#[repr(C)]
pub(super) struct MitoEngineAnalyzeOptionsV1 {
    pub struct_size: u32,
    pub abi_version: u32,
    pub filter_numt: bool,
    pub threads: usize,
    pub min_mapping_quality: u8,
    pub min_base_quality: u8,
    pub excluded_snp_flags: u16,
    pub numt_threshold: f64,
    pub allow_development_tags: bool,
    pub emit_evidence_graph: bool,
    pub max_evidence_observations: usize,
    pub max_phase_links: usize,
    pub max_phase_work: usize,
    pub max_phase_molecule_references: usize,
    pub max_result_bytes: usize,
    pub evidence_page_size: usize,
    pub molecule_id_tag: *const c_char,
    pub umi_tag: *const c_char,
    pub duplex_tag: *const c_char,
    pub min_architecture_molecules: usize,
    pub max_candidate_architectures: usize,
    pub architecture_max_distance: f64,
    pub architecture_ambiguity_margin: f64,
    pub architecture_min_overlap_fraction: f64,
    pub architecture_consensus_fraction: f64,
    pub architecture_optional_fraction: f64,
    pub architecture_stability_replicates: usize,
    pub architecture_seed: u64,
    pub should_cancel: Option<unsafe extern "C" fn(*mut c_void) -> bool>,
    pub cancel_user_data: *mut c_void,
}

unsafe extern "C" {
    pub(super) fn mito_engine_new() -> *mut c_void;
    pub(super) fn mito_engine_delete(engine: *mut c_void);
    pub(super) fn mito_engine_has_htslib() -> bool;
    pub(super) fn mito_engine_version() -> *const c_char;
    pub(super) fn mito_engine_schema_version() -> *const c_char;
    pub(super) fn mito_engine_error_schema_version() -> *const c_char;
    pub(super) fn mito_engine_analyze_with_options_v1(
        engine: *mut c_void,
        input_path: *const c_char,
        ref_path: *const c_char,
        options: *const MitoEngineAnalyzeOptionsV1,
    ) -> *const c_char;
    pub(super) fn mito_engine_get_last_error() -> *const c_char;
    pub(super) fn mito_engine_get_last_error_code() -> *const c_char;
    pub(super) fn mito_engine_free_string(value: *const c_char);
}
