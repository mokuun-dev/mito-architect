use crate::args::{AnalyzeArgs, UpdateClinicalArgs};
use crate::clinical::update_clinical;
use crate::exports::{
    build_provenance_manifest, export_evidence_pages, render_variant_tsv, render_vcf,
    ExportArtifact,
};
use crate::report;
use crate::storage::{write_bgzip_vcf_with_tabix, write_bytes_atomic, write_json_atomic};
use anyhow::Context;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use mito_ffi::AnalyzeOptions;
use mito_ffi::MitoEngine;
use std::time::Duration;

pub(super) fn analyze(args: AnalyzeArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("input file does not exist: {}", args.input.display());
    }
    if let Some(reference) = &args.reference {
        if !reference.exists() {
            anyhow::bail!("reference file does not exist: {}", reference.display());
        }
    }
    if args.update_clinical {
        let cache_path = update_clinical(UpdateClinicalArgs {
            source: None,
            output: None,
            clinvar_live: false,
            clinvar_gz: None,
        })?;
        // SAFETY: the CLI is single-threaded before the C++ engine is created.
        unsafe {
            std::env::set_var("MITO_CLINICAL_ANNOTATIONS", &cache_path);
        }
    }

    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::with_template("{spinner:.cyan} {msg}")
            .unwrap_or_else(|_| ProgressStyle::default_spinner()),
    );
    progress.enable_steady_tick(Duration::from_millis(80));
    progress.set_message("initializing C++ analysis engine");

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    progress.set_message("analyzing reads");
    let emit_evidence_graph = args.emit_evidence_graph
        || args.tsv.is_some()
        || args.bgzip_vcf.is_some()
        || args.evidence_pages_dir.is_some();
    let result = engine
        .analyze_with_options(
            &args.input,
            args.reference.as_deref(),
            AnalyzeOptions {
                filter_numt: args.filter_numt,
                threads: args.threads,
                min_mapping_quality: args.min_mapping_quality,
                min_base_quality: args.min_base_quality,
                excluded_snp_flags: args.excluded_snp_flags,
                numt_threshold: args.numt_threshold,
                allow_development_tags: args.allow_development_tags,
                emit_evidence_graph,
                max_evidence_observations: args.max_evidence_observations,
                max_phase_links: args.max_phase_links,
                max_phase_work: args.max_phase_work,
                max_phase_molecule_references: args.max_phase_molecule_references,
                max_result_bytes: args.max_result_bytes,
                evidence_page_size: args.evidence_page_size,
                min_architecture_molecules: args.min_architecture_molecules,
                max_candidate_architectures: args.max_candidate_architectures,
                architecture_max_distance: args.architecture_max_distance,
                architecture_ambiguity_margin: args.architecture_ambiguity_margin,
                architecture_min_overlap_fraction: args.architecture_min_overlap_fraction,
                architecture_consensus_fraction: args.architecture_consensus_fraction,
                architecture_optional_fraction: args.architecture_optional_fraction,
                architecture_stability_replicates: args.architecture_stability_replicates,
                architecture_seed: args.architecture_seed,
                molecule_id_tag: args.molecule_id_tag.clone(),
                umi_tag: args.umi_tag.clone(),
                duplex_tag: args.duplex_tag.clone(),
            },
        )
        .context("analysis failed")?;

    let mut exported = Vec::<ExportArtifact>::new();

    if let Some(directory) = &args.evidence_pages_dir {
        progress.set_message("writing evidence page resources");
        export_evidence_pages(&result, directory)?;
        let manifest = directory.join("manifest.json");
        exported.push(ExportArtifact::from_path("evidence_manifest", &manifest)?);
    }

    let vcf = if args.vcf.is_some() || args.bgzip_vcf.is_some() {
        progress.set_message("writing VCF export");
        Some(render_vcf(&result)?)
    } else {
        None
    };
    if let (Some(vcf_path), Some(vcf)) = (&args.vcf, vcf.as_deref()) {
        write_bytes_atomic(vcf_path, vcf.as_bytes())?;
        exported.push(ExportArtifact::from_path("vcf", vcf_path)?);
    }
    if let (Some(vcf_path), Some(vcf)) = (&args.bgzip_vcf, vcf.as_deref()) {
        progress.set_message("writing bgzip VCF and tabix index");
        let index_path = write_bgzip_vcf_with_tabix(vcf_path, vcf)?;
        exported.push(ExportArtifact::from_path("vcf_bgzip", vcf_path)?);
        exported.push(ExportArtifact::from_path("vcf_tabix", &index_path)?);
    }

    if let Some(tsv_path) = &args.tsv {
        progress.set_message("writing unified variant TSV");
        let tsv = render_variant_tsv(&result)?;
        write_bytes_atomic(tsv_path, tsv.as_bytes())?;
        exported.push(ExportArtifact::from_path("variant_tsv", tsv_path)?);
    }

    if !args.json {
        progress.set_message("rendering standalone report");
        let html = report::render_report(&result)?;
        write_bytes_atomic(&args.output, html.as_bytes())?;
        exported.push(ExportArtifact::from_path("html", &args.output)?);
    }

    if let Some(manifest_path) = &args.provenance_manifest {
        progress.set_message("writing provenance manifest");
        let manifest = build_provenance_manifest(&result, &args, &exported)?;
        write_json_atomic(manifest_path, &manifest)?;
    }

    if args.json {
        progress.finish_and_clear();
        println!("{result}");
    } else {
        progress.finish_with_message(format!("wrote {}", args.output.display()));
    }

    Ok(())
}
