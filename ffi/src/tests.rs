use super::{AnalyzeOptions, MitoEngine};
use std::fs;
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};

#[test]
fn cancelled_analysis_returns_error() {
    let input_path = std::env::temp_dir().join(format!(
        "mito_ffi_cancel_{}_{}.fastq",
        std::process::id(),
        unique_suffix()
    ));
    let mut input = fs::File::create(&input_path).expect("temp input");
    writeln!(input, "@read-1").expect("write header");
    writeln!(input, "GATCACAGGT").expect("write sequence");
    writeln!(input, "+").expect("write plus");
    writeln!(input, "IIIIIIIIII").expect("write quality");
    drop(input);

    let cancel = AtomicBool::new(true);
    let engine = MitoEngine::new().expect("engine");
    let error = engine
        .analyze_with_cancel_flag(
            &input_path,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads: 1,
                ..AnalyzeOptions::default()
            },
            &cancel,
        )
        .expect_err("analysis should be cancelled");
    assert_eq!(error.code, "MITO-E1501");
    assert_eq!(error.message, "analysis cancelled");
    assert!(cancel.load(Ordering::Relaxed));
    let _ = fs::remove_file(input_path);
}

#[test]
fn evidence_graph_is_opt_in_and_bounded() {
    let input_path = temp_path("evidence-graph", "sam");
    fs::write(
            &input_path,
            "@SQ\tSN:NC_012920.1\tLN:16569\nalt\t0\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\tIII\nreference\t0\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAT\tIII\n",
        )
        .expect("evidence input");

    let engine = MitoEngine::new().expect("engine");
    let json = engine
        .analyze_with_options(
            &input_path,
            None,
            AnalyzeOptions {
                emit_evidence_graph: true,
                evidence_page_size: 1,
                ..AnalyzeOptions::default()
            },
        )
        .expect("schema 0.6 result");
    assert!(json.contains("\"schema_version\":\"0.6\""));
    assert!(json.contains("\"missing_pair_state\":\"NOT_CALLABLE\""));
    assert!(json.contains("\"layout\":\"paged_columnar_molecule_event\""));
    assert!(json.contains("\"observation_page_size\":1"));

    let error = engine
        .analyze_with_options(
            &input_path,
            None,
            AnalyzeOptions {
                emit_evidence_graph: true,
                max_evidence_observations: 1,
                ..AnalyzeOptions::default()
            },
        )
        .expect_err("observation cap must fail closed");
    assert_eq!(error.code, "MITO-E1601");
    let _ = fs::remove_file(input_path);
}

#[test]
fn protocol_tags_cross_the_versioned_ffi_boundary() {
    let input_path = temp_path("protocol-tags", "sam");
    fs::write(
            &input_path,
            "@SQ\tSN:NC_012920.1\tLN:16569\nread-a\t0\tNC_012920.1\t1\t60\t3M\t*\t0\t0\tGAA\tIII\tMI:Z:M1\tRX:Z:AAA\tDX:Z:duplex\n",
        )
        .expect("protocol input");
    let engine = MitoEngine::new().expect("engine");
    let json = engine
        .analyze_with_options(
            &input_path,
            None,
            AnalyzeOptions {
                emit_evidence_graph: true,
                molecule_id_tag: "MI".to_owned(),
                umi_tag: "RX".to_owned(),
                duplex_tag: "DX".to_owned(),
                ..AnalyzeOptions::default()
            },
        )
        .expect("tagged schema 0.6 result");
    assert!(json.contains("\"identity_policy\":\"sam_tag:MI\""));
    assert!(json.contains("\"molecule_id_value\":\"M1\""));
    assert!(json.contains("\"UMI_RECORDED\""));
    assert!(json.contains("\"DUPLEX_METADATA_RECORDED\""));
    let _ = fs::remove_file(input_path);
}

#[test]
fn external_failures_preserve_stable_error_codes() {
    let valid_input = temp_path("valid", "fastq");
    fs::write(&valid_input, "@read-1\nGATCACAGGT\n+\nIIIIIIIIII\n").expect("valid input");
    let malformed_input = temp_path("malformed", "fastq");
    fs::write(&malformed_input, "@read-1\nGATC\nnot-plus\nIIII\n").expect("malformed input");
    let fasta_input = temp_path("reads", "fasta");
    fs::write(&fasta_input, ">read-1\nGATC\n").expect("FASTA input");
    let invalid_reference = temp_path("invalid-reference", "fasta");
    fs::write(&invalid_reference, "GATC\n").expect("invalid reference");
    let empty_input = temp_path("empty", "fastq");
    fs::write(&empty_input, "").expect("empty input");
    let missing_input = temp_path("missing", "fastq");
    let missing_reference = temp_path("missing-reference", "fasta");

    let engine = MitoEngine::new().expect("engine");
    let cases = [
        (
            &malformed_input,
            None,
            AnalyzeOptions::default(),
            "MITO-E1103",
        ),
        (&fasta_input, None, AnalyzeOptions::default(), "MITO-E1102"),
        (&empty_input, None, AnalyzeOptions::default(), "MITO-E1104"),
        (
            &missing_input,
            None,
            AnalyzeOptions::default(),
            "MITO-E1101",
        ),
        (
            &valid_input,
            Some(invalid_reference.as_path()),
            AnalyzeOptions::default(),
            "MITO-E1202",
        ),
        (
            &valid_input,
            Some(missing_reference.as_path()),
            AnalyzeOptions::default(),
            "MITO-E1201",
        ),
        (
            &valid_input,
            None,
            AnalyzeOptions {
                numt_threshold: f64::NAN,
                ..AnalyzeOptions::default()
            },
            "MITO-E1001",
        ),
    ];
    for (input, reference, options, expected_code) in cases {
        let error = engine
            .analyze_with_options(input, reference, options)
            .expect_err("negative case must fail");
        assert_eq!(error.code, expected_code, "unexpected error: {error}");
    }

    for path in [
        valid_input,
        malformed_input,
        fasta_input,
        invalid_reference,
        empty_input,
    ] {
        let _ = fs::remove_file(path);
    }
}

fn temp_path(label: &str, extension: &str) -> std::path::PathBuf {
    std::env::temp_dir().join(format!(
        "mito_ffi_{label}_{}_{}.{}",
        std::process::id(),
        unique_suffix(),
        extension
    ))
}

fn unique_suffix() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or_default()
}

#[test]
fn nul_bytes_are_rejected_before_native_calls() {
    use std::path::Path;
    let engine = MitoEngine::new().unwrap();
    let error = engine.analyze(Path::new("bad\0path"), None).unwrap_err();
    assert_eq!(error.code, "MITO-E1001");
    let error = engine
        .analyze_with_options(
            Path::new("unused.fastq"),
            None,
            AnalyzeOptions {
                molecule_id_tag: "M\0".into(),
                ..Default::default()
            },
        )
        .unwrap_err();
    assert_eq!(error.code, "MITO-E1001");
}

#[test]
fn engine_can_move_threads_and_recover_after_cancellation() {
    let engine = MitoEngine::new().unwrap();
    std::thread::spawn(move || {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../fixtures/tiny.fastq");
        let flag = AtomicBool::new(true);
        assert_eq!(
            engine
                .analyze_with_cancel_flag(&path, None, Default::default(), &flag)
                .unwrap_err()
                .code,
            "MITO-E1501"
        );
        flag.store(false, Ordering::Relaxed);
        let json = engine
            .analyze_with_cancel_flag(&path, None, Default::default(), &flag)
            .unwrap();
        assert!(json.contains("\"input_molecules\":4"));
    })
    .join()
    .unwrap();
}

#[cfg(unix)]
#[test]
fn native_output_escapes_invalid_input_bytes_as_valid_json() {
    let path = temp_path("invalid-utf8", "fastq");
    fs::write(&path, b"@read-\xff\nGATC\n+\nIIII\n").unwrap();
    let result = MitoEngine::new().unwrap().analyze(&path, None);
    fs::remove_file(path).unwrap();
    let json = result.expect("invalid input bytes must not corrupt native JSON");
    assert!(json.contains("read-\\u00ff"));
}

#[cfg(target_pointer_width = "64")]
#[test]
fn options_struct_matches_the_c_layout() {
    use crate::raw::MitoEngineAnalyzeOptionsV1;
    use std::mem::offset_of;
    // Mirrors check_options_layout in ffi/tests/c_api_contract.c (LP64).
    let offsets = [
        (
            "struct_size",
            offset_of!(MitoEngineAnalyzeOptionsV1, struct_size),
            0,
        ),
        (
            "abi_version",
            offset_of!(MitoEngineAnalyzeOptionsV1, abi_version),
            4,
        ),
        (
            "filter_numt",
            offset_of!(MitoEngineAnalyzeOptionsV1, filter_numt),
            8,
        ),
        (
            "threads",
            offset_of!(MitoEngineAnalyzeOptionsV1, threads),
            16,
        ),
        (
            "min_mapping_quality",
            offset_of!(MitoEngineAnalyzeOptionsV1, min_mapping_quality),
            24,
        ),
        (
            "min_base_quality",
            offset_of!(MitoEngineAnalyzeOptionsV1, min_base_quality),
            25,
        ),
        (
            "excluded_snp_flags",
            offset_of!(MitoEngineAnalyzeOptionsV1, excluded_snp_flags),
            26,
        ),
        (
            "numt_threshold",
            offset_of!(MitoEngineAnalyzeOptionsV1, numt_threshold),
            32,
        ),
        (
            "allow_development_tags",
            offset_of!(MitoEngineAnalyzeOptionsV1, allow_development_tags),
            40,
        ),
        (
            "emit_evidence_graph",
            offset_of!(MitoEngineAnalyzeOptionsV1, emit_evidence_graph),
            41,
        ),
        (
            "max_evidence_observations",
            offset_of!(MitoEngineAnalyzeOptionsV1, max_evidence_observations),
            48,
        ),
        (
            "max_phase_links",
            offset_of!(MitoEngineAnalyzeOptionsV1, max_phase_links),
            56,
        ),
        (
            "max_phase_work",
            offset_of!(MitoEngineAnalyzeOptionsV1, max_phase_work),
            64,
        ),
        (
            "max_phase_molecule_references",
            offset_of!(MitoEngineAnalyzeOptionsV1, max_phase_molecule_references),
            72,
        ),
        (
            "max_result_bytes",
            offset_of!(MitoEngineAnalyzeOptionsV1, max_result_bytes),
            80,
        ),
        (
            "evidence_page_size",
            offset_of!(MitoEngineAnalyzeOptionsV1, evidence_page_size),
            88,
        ),
        (
            "molecule_id_tag",
            offset_of!(MitoEngineAnalyzeOptionsV1, molecule_id_tag),
            96,
        ),
        (
            "umi_tag",
            offset_of!(MitoEngineAnalyzeOptionsV1, umi_tag),
            104,
        ),
        (
            "duplex_tag",
            offset_of!(MitoEngineAnalyzeOptionsV1, duplex_tag),
            112,
        ),
        (
            "min_architecture_molecules",
            offset_of!(MitoEngineAnalyzeOptionsV1, min_architecture_molecules),
            120,
        ),
        (
            "max_candidate_architectures",
            offset_of!(MitoEngineAnalyzeOptionsV1, max_candidate_architectures),
            128,
        ),
        (
            "architecture_max_distance",
            offset_of!(MitoEngineAnalyzeOptionsV1, architecture_max_distance),
            136,
        ),
        (
            "architecture_ambiguity_margin",
            offset_of!(MitoEngineAnalyzeOptionsV1, architecture_ambiguity_margin),
            144,
        ),
        (
            "architecture_min_overlap_fraction",
            offset_of!(
                MitoEngineAnalyzeOptionsV1,
                architecture_min_overlap_fraction
            ),
            152,
        ),
        (
            "architecture_consensus_fraction",
            offset_of!(MitoEngineAnalyzeOptionsV1, architecture_consensus_fraction),
            160,
        ),
        (
            "architecture_optional_fraction",
            offset_of!(MitoEngineAnalyzeOptionsV1, architecture_optional_fraction),
            168,
        ),
        (
            "architecture_stability_replicates",
            offset_of!(
                MitoEngineAnalyzeOptionsV1,
                architecture_stability_replicates
            ),
            176,
        ),
        (
            "architecture_seed",
            offset_of!(MitoEngineAnalyzeOptionsV1, architecture_seed),
            184,
        ),
        (
            "should_cancel",
            offset_of!(MitoEngineAnalyzeOptionsV1, should_cancel),
            192,
        ),
        (
            "cancel_user_data",
            offset_of!(MitoEngineAnalyzeOptionsV1, cancel_user_data),
            200,
        ),
        (
            "max_input_records",
            offset_of!(MitoEngineAnalyzeOptionsV1, max_input_records),
            208,
        ),
        (
            "max_input_bases",
            offset_of!(MitoEngineAnalyzeOptionsV1, max_input_bases),
            216,
        ),
        (
            "allow_external_cram_references",
            offset_of!(MitoEngineAnalyzeOptionsV1, allow_external_cram_references),
            224,
        ),
    ];
    for (field, actual, expected) in offsets {
        assert_eq!(actual, expected, "offset of {field}");
    }
    assert_eq!(std::mem::size_of::<MitoEngineAnalyzeOptionsV1>(), 232);
}
