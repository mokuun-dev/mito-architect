use crate::architecture_contract::validate_architecture_contract;
use crate::exports::{render_variant_tsv, render_vcf};
use crate::fixtures::{assert_aux_present, assert_mapq_present, assert_sv_present};
use crate::projection::materialize_observations;

#[test]
fn renders_vcf_with_locus_callable_heteroplasmy() {
    let json = r#"{
      "metadata": {
        "schema_version": "0.5",
        "sv_event_schema_version": "1.0",
        "complex_sv_event_schema_version": "1.0",
        "clinical_annotation_schema_version": "1.0",
        "engine_version": "0.5.0-dev",
        "reference_length": 16569,
        "sample": "sample A",
        "reference_accession": "NC_012920.1",
        "resources": []
      },
      "filter_stats": {
        "passed_reads": 40,
        "numt_assessment": { "mode": "competitive_alignment" }
      },
      "variants": [
        {
          "position": 3243,
          "ref": "A",
          "alt": "G",
          "alt_depth": 2,
          "callable_depth": 3,
          "ci95_low": 0.207660,
          "ci95_high": 0.938508,
          "gene": "MT-TL1",
          "annotation": { "pathogenicity": "pathogenic" }
        }
      ],
      "reads": [],
      "svs": [],
      "complex_events": [],
      "clusters": []
    }"#;

    let vcf = render_vcf(json).expect("VCF should render");
    assert!(vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_A"));
    assert!(vcf.contains("NC_012920.1\t3243\t.\tA\tG\t.\tPASS\tAC=2;DP=3;HF=0.666667;GENE=MT-TL1;CLNSIG=pathogenic;HF_CI95=0.207660,0.938508\tGT:HF\t0/1:0.666667"));
}

#[test]
fn renders_schema_0_6_small_indel_with_stable_event_id() {
    let json = r#"{
      "metadata": {
        "schema_version": "0.6",
        "sv_event_schema_version": "1.0",
        "complex_sv_event_schema_version": "1.0",
        "clinical_annotation_schema_version": "1.0",
        "engine_version": "0.5.0-dev",
        "reference_length": 16569,
        "sample": "phase",
        "reference_accession": "NC_012920.1",
        "resources": []
      },
      "filter_stats": {
        "passed_reads": 2,
        "numt_assessment": { "mode": "mt_only_or_unknown" }
      },
      "evidence_encoding": {
        "layout": "paged_columnar_molecule_event",
        "observation_storage": "embedded_columnar_pages",
        "missing_pair_state": "NOT_CALLABLE",
        "phase_molecule_policy": "evidence_eligible_only",
        "phase_molecule_reference": "molecules[].index",
        "phase_null_model": "independent_marginals_within_jointly_callable",
        "observation_limit": 10,
        "observation_count": 2,
        "observation_page_size": 2,
        "observation_page_count": 1,
        "phase_link_limit": 10
      },
      "alignments": [
        {"id":"alignment:0","molecule_id":"m1"},
        {"id":"alignment:1","molecule_id":"m2"}
      ],
      "molecules": [
        {"id":"m1","index":0,"identity_policy":"sam_qname","source_qnames":["m1"],"protocol_metadata":{},"analysis_eligible":true,"evidence_eligible":true,"exclusion_reasons":[],"representative_alignment_id":"alignment:0","alignment_ids":["alignment:0"],"query_length":12,"mean_base_quality":40,"mapping_quality":60,"numt_score":0,"numt_evidence":[],"architecture_assignment":{"status":"UNASSIGNED","architecture_id":null,"candidate_architecture_ids":[]},"alternate_event_ids":["indel:insertion:5:G"],"evidence_state_counts":{"alternate":1,"reference":0,"event_absent":0,"low_quality":0,"conflict":0}},
        {"id":"m2","index":1,"identity_policy":"sam_qname","source_qnames":["m2"],"protocol_metadata":{},"analysis_eligible":true,"evidence_eligible":true,"exclusion_reasons":[],"representative_alignment_id":"alignment:1","alignment_ids":["alignment:1"],"query_length":12,"mean_base_quality":40,"mapping_quality":60,"numt_score":0,"numt_evidence":[],"architecture_assignment":{"status":"UNASSIGNED","architecture_id":null,"candidate_architecture_ids":[]},"alternate_event_ids":[],"evidence_state_counts":{"alternate":0,"reference":0,"event_absent":1,"low_quality":0,"conflict":0}}
      ],
      "callability": [
        {"molecule_id":"m1","ranges":[{"start":1,"end":12}],"alignments":[{"alignment_id":"alignment:0"}]},
        {"molecule_id":"m2","ranges":[{"start":1,"end":12}],"alignments":[{"alignment_id":"alignment:1"}]}
      ],
      "events": [{
        "id":"indel:insertion:5:G",
        "index":0,
        "type":"SMALL_INSERTION",
        "start":5,
        "end":5,
        "ref":"A",
        "alt":"AG",
        "negative_evidence_rule":"same_fragment_callable_reference_adjacency",
        "assessability":"REFERENCE_AND_ALTERNATE",
        "supporting_molecule_ids":["m1"],
        "evidence_counts":{"alternate":1,"reference":0,"event_absent":1,"callable":2,"low_quality":0,"conflict":0}
      }],
      "observation_pages": [
        {
          "index": 0,
          "offset": 0,
          "count": 2,
          "columns": {
            "molecule_id": ["m1", "m2"],
            "event_id": ["indel:insertion:5:G", "indel:insertion:5:G"],
            "alignment_id": ["alignment:0", "alignment:1"],
            "state": ["ALTERNATE", "EVENT_ABSENT"],
            "observed_allele": ["G", null],
            "base_quality": [40, null],
            "mapping_quality": [60, 60],
            "strand": ["+", "+"],
            "evidence_source": ["cigar_small_indel", "callable_reference_path"],
            "read_position": [0, null]
          }
        }
      ],
      "phase_links": [],
      "architecture_inference": {
        "status":"NO_SUPPORTED_SEEDS",
        "eligible_molecules":2,
        "assigned_molecules":0,
        "ambiguous_molecules":0,
        "unassigned_molecules":2,
        "candidate_count":0
      },
      "architectures": [],
      "variants": [{
        "event_id":"indel:insertion:5:G",
        "type":"SMALL_INSERTION",
        "position":5,
        "start":5,
        "end":5,
        "ref":"A",
        "alt":"AG",
        "normalization":"rcrs_3prime_small_indel_v1",
        "vcf_position":5,
        "vcf_representable":true,
        "alt_depth":1,
        "ref_depth":1,
        "other_depth":0,
        "event_absent_depth":1,
        "low_quality_depth":0,
        "conflict_depth":0,
        "callable_depth":2,
        "heteroplasmy":0.5,
        "ci95_low":0.094531,
        "ci95_high":0.905469,
        "filter_status":"NOT_CALIBRATED",
        "qc_flags":["NUMT_NOT_ASSESSABLE","SINGLE_STRAND_ALT_SUPPORT"],
        "numt_assessability":"NOT_ASSESSABLE",
        "supporting_molecule_ids":["m1"],
        "mapping_quality":{"alternate":{"mean":60}},
        "allele_quality":{"alternate":{"mean_phred":40}},
        "strand_support":{"alt_forward":1,"alt_reverse":0,"ref_forward":1,"ref_reverse":0}
      }],
      "reads": [],
      "svs": [],
      "complex_events": [],
      "clusters": []
    }"#;

    let vcf = render_vcf(json).expect("schema 0.6 VCF should render");
    assert!(vcf.contains("##INFO=<ID=EVENT_TYPE"));
    assert!(vcf.contains(
        "NC_012920.1\t5\tindel:insertion:5:G\tA\tAG\t.\t.\tAC=1;DP=2;HF=0.500000;AD=1,1;ODC=0;LOWQ=0;CONFLICT=0;MOLECULE_SUPPORT=1;MQ=60.0000;BQ=40.0000;STRAND_SUPPORT=1,0,1,0;NUMT_ASSESSABLE=NOT_ASSESSABLE;QC_FLAGS=NUMT_NOT_ASSESSABLE,SINGLE_STRAND_ALT_SUPPORT;NORMALIZATION=rcrs_3prime_small_indel_v1;HF_CI95=0.094531,0.905469;EVENT_TYPE=SMALL_INSERTION\tGT:DP:AD:HF\t.:2:1,1:0.500000"
    ));
    let tsv = render_variant_tsv(json).expect("schema 0.6 TSV should render");
    assert!(tsv
        .contains("indel:insertion:5:G\tSMALL_INSERTION\t5\t5\t5\tA\tAG\t1\t1\t0\t2\t0.50000000"));
}

#[test]
fn materializes_columnar_observation_pages_and_rejects_cardinality_drift() {
    let mut value = serde_json::json!({
        "evidence_encoding": {
            "observation_count": 2,
            "observation_page_size": 2,
            "observation_page_count": 1
        },
        "observation_pages": [{
            "index": 0,
            "offset": 0,
            "count": 2,
            "columns": {
                "molecule_id": ["m1", "m2"],
                "event_id": ["e1", "e2"],
                "alignment_id": ["alignment:0", "alignment:1"],
                "state": ["ALTERNATE", "REFERENCE"],
                "observed_allele": ["G", "A"],
                "base_quality": [40, 39],
                "mapping_quality": [60, 50],
                "strand": ["+", "-"],
                "evidence_source": ["aligned_base", "aligned_base"],
                "read_position": [0, 1]
            }
        }]
    });
    let rows = materialize_observations(&value).expect("page should materialize");
    assert_eq!(rows.len(), 2);
    assert_eq!(
        rows[1].get("id").and_then(serde_json::Value::as_str),
        Some("observation:1")
    );

    value["observation_pages"][0]["columns"]["strand"] = serde_json::json!(["+"]);
    let error = materialize_observations(&value).expect_err("column drift must fail closed");
    assert!(error
        .to_string()
        .contains("column strand has the wrong length"));
}

#[test]
fn validates_fixture_fields_beyond_snp() {
    let json: serde_json::Value = serde_json::from_str(
        r#"{
          "coverage_metrics": {
            "mapping_quality_histogram": [
              { "mapq": 42, "count": 1 }
            ]
          },
          "svs": [
            { "id": "deletion:110-121" }
          ],
          "reads": [
            { "aux_tags": { "NM": "1", "MD": "2T7" } }
          ]
        }"#,
    )
    .expect("fixture JSON should parse");

    assert_sv_present(&json, "deletion:110-121").expect("SV should validate");
    assert_mapq_present(&json, 42).expect("MAPQ should validate");
    assert_aux_present(&json, "NM=1").expect("aux tag should validate");
    assert_aux_present(&json, "MD=2T7").expect("aux tag should validate");
}

#[test]
fn architecture_event_profiles_fail_closed_on_count_drift() {
    let mut result: serde_json::Value = serde_json::from_str(include_str!(
        "../../fixtures/truth_architectures.expected.json"
    ))
    .expect("architecture golden should parse");
    let event_ids = result["events"]
        .as_array()
        .expect("events should be an array")
        .iter()
        .map(|event| {
            event["id"]
                .as_str()
                .expect("event ID should be a string")
                .to_owned()
        })
        .collect::<std::collections::BTreeSet<_>>();
    let eligible_molecules = result["molecules"]
        .as_array()
        .expect("molecules should be an array")
        .iter()
        .filter(|molecule| molecule["evidence_eligible"].as_bool() == Some(true))
        .map(|molecule| {
            (
                molecule["index"].as_u64().expect("index should be numeric") as usize,
                molecule["id"]
                    .as_str()
                    .expect("molecule ID should be a string")
                    .to_owned(),
            )
        })
        .collect::<std::collections::BTreeMap<_, _>>();
    validate_architecture_contract(&result, &event_ids, &eligible_molecules)
        .expect("architecture golden should validate");

    result["architectures"][0]["event_profile"][0]["not_callable"] = serde_json::json!(1);
    let error = validate_architecture_contract(&result, &event_ids, &eligible_molecules)
        .expect_err("architecture member-count drift must fail closed");
    assert!(error.to_string().contains("does not partition members"));
}
