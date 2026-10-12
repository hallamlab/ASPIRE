#!/usr/bin/env python3
"""Generate the exhaustive ASPIRE YAML parameter catalogue."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterator

import yaml


ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "asv_pipeline_nextflow.yml"
OUTPUT = ROOT / "docs" / "CONFIG_PARAMETERS.md"

EXACT = {
    "optional.power_analysis.sample_sizes_cancer": "Cancer-patient counts to simulate. Controls retain the observed pilot count until cancer count exceeds it, then use equal counts per group.",
    "optional.power_analysis.sample_sizes_stype": "Patient-profile counts to simulate for paired sample-type comparisons; not library counts or sequencing depth.",
    "optional.power_analysis.workers": "Maximum simulation worker processes. Null inherits core.resources.threads; explicit values are capped at that budget and available CPUs. One runs serially.",
    "optional.power_analysis.n_simulations": "Number of repeated simulations per scenario and patient-count setting; controls Monte Carlo precision.",
    "optional.power_analysis.n_perm": "Number of permutations within each permutation test; separate from the simulation count.",
    "optional.power_analysis.alpha": "Significance cutoff used to count detections in power simulations.",
    "standard.mito.mitomaster_failure_policy": "After bounded API retries, fail stops the process; continue retains successful responses and runs local BLAST even if all API requests fail. Publishes status JSON and failed-chunk TSV; unavailable API evidence is not a completed negative screen.",
    "optional.sankey.arrangement": "Vertical node interaction: snap returns nodes to their ordered initial positions on release; freeform and perpendicular retain constrained vertical positions; fixed disables dragging. Retained/removed lanes and within-column order are preserved in all modes.",
    "optional.analysis_cohort.sample_col": "Sample identifier column used to synchronize downstream metadata, counts and long tables.",
    "optional.analysis_cohort.group_col": "Metadata column whose exact values define study groups to exclude from downstream analysis.",
    "optional.analysis_cohort.exclude_groups": "List of study groups retained in QC, decontamination, metadata plots, batch preparation and diversity, but excluded from every other sample-based optional analysis. Empty disables cohort selection. Does not change control roles or upstream abundance/prevalence denominators.",
    "optional.voc_correlation.isa_exclude_nondistinct": "Exclude ISA memberships containing every isa_all_type_groups value from all ISA-specific VOC plots and tables. Excluded memberships are audited; general ASV correlation tables remain available.",
    "core.primer_trimming.screen_error_rate": "Maximum substitution rate for complete primers during screening (0 to 0.25), without indels; 0 requires exact IUPAC matches. Independent of Cutadapt error_rate. Allowed substitutions are floor(rate times primer length).",
    "environments.primer_trimming": "Conda environment YAML for paired Cutadapt primer removal and cohort checks.",
    'core.primer_trimming.enabled': 'Run paired Cutadapt primer removal and cohort family validation before fastp. Requires paired reads and all four fastp fixed trimming values set to zero.',
    'core.primer_trimming.primers': 'Candidate primer families, each with unique name, forward and reverse IUPAC DNA sequences. Restrict to the known assay when possible; candidate screening does not discover unknown primers.',
    'core.primer_trimming.sample_reads': 'Number of initial read pairs screened per sample; positive integer. Screening is sequential, not random.',
    'core.primer_trimming.max_prefix': 'Maximum number of bases allowed before a 5-prime primer; integer from 0 to 100. The prefix is removed with the primer.',
    'core.primer_trimming.min_pair_fraction': 'Minimum fraction of screened pairs supporting a candidate primer pair (greater than 0 and at most 1). Insufficient support stops the sample.',
    'core.primer_trimming.min_family_fraction': 'Minimum fraction of screened pairs supporting a family/orientation for selection (greater than 0 and at most 1).',
    'core.primer_trimming.min_family_reads': 'Minimum supporting screened pairs for each selected family/orientation; positive integer.',
    'core.primer_trimming.error_rate': 'Maximum Cutadapt primer mismatch rate (0 to 0.25); indels are disabled. Screening has a separate screen_error_rate setting.',
    'core.primer_trimming.end_overlap': 'Minimum overlap for optional opposite-primer read-through removal at the 3-prime end; positive integer.',
    'core.primer_trimming.minimum_length': 'Minimum length of each primer-trimmed mate; positive integer. Pairs with a shorter mate enter the too-short audit bin.',

    "standard.filter_counts.min_prevalence_fraction": "Minimum fraction of retained biological samples with nonzero counts, independent of RA. General and mock default 0.05 (5%); 0 disables it. Positive values require control decontamination. Includes zero-depth columns remaining after feature removal.",
    "optional.taxonomy_patient_aware.min_prevalence_fraction": "Fraction of nonzero patient profiles retained for taxonomic tests; for sample-type contrasts use pooled patient-by-type profiles. 0.1 means 10%.",
    "optional.measurement_association.min_prevalence_fraction": "Fraction of matched measurement/count samples with nonzero ASV counts. 0 disables this filter; this is not the final biological-table filter.",
    "standard.metadata_plots.input_table": "Must be filtered: PLOT_METADATA receives the complete FILTER_ASVS result; pre-filter microbial counts are audit-only.",
    'core.control_decontam.enabled': 'Run independent TECH/BIO prevalence tests before biological feature filtering.',
    'core.control_decontam.class_col': 'Metadata column defining sample class. Filenames and sample IDs are never used to infer class.',
    'core.control_decontam.min_biological_reads': 'Minimum post-QC count sum for biological samples only. Biological/sample and technical controls are exempt.',
    'core.control_decontam.biological_labels': 'Metadata values identifying biological samples.',
    'core.control_decontam.technical_labels': 'Metadata values identifying water, extraction and library negatives.',
    'core.control_decontam.bio_control_labels': 'Metadata values identifying biological/sample controls, such as scope flushes.',
    'core.control_decontam.positive_labels': 'Metadata values identifying positive controls excluded from both tests.',
    'core.control_decontam.technical_enabled': 'Run technical-negative versus QC-passing biological prevalence test.',
    'core.control_decontam.bio_control_enabled': 'Run biological-control versus QC-passing biological prevalence test.',
    'core.control_decontam.technical_score_threshold': 'TECH decontam score cutoff, strictly between 0 and 1; flag scores below it. Default 0.1 is not 0.1% or minimum sample prevalence.',
    'core.control_decontam.bio_control_score_threshold': 'BIO decontam score cutoff, strictly between 0 and 1; flag scores below it. Independent of TECH and of abundance/prevalence filters.',
    "optional.batch_correction.correction_policy": "Select downstream counts: always uses corrected counts, never uses raw counts, and auto applies the count-space/batch/biological preservation gates. The template preserves the original always-correct behavior.",
    "optional.batch_correction.auto_min_sample_rho": "Minimum median per-sample Spearman count-preservation correlation required by automatic correction selection.",
    "optional.batch_correction.auto_min_bray_rho": "Minimum Bray-Curtis distance correlation required by automatic correction selection.",
    "optional.spieceasi.pulsar_criterion": "Model selection: stars uses ordinary StARS; bstars requests bounded StARS via lower/upper bounds. The template preserves ordinary StARS.",
    "optional.measurement_association.measurement_table": "Optional sample-level measurement table; when absent, eligible numeric metadata columns supply measurements.",
    "optional.measurement_association.metadata_join_cols": "Metadata join keys paired positionally with measurement_join_cols; empty lists use the configured sample identifiers.",
    "optional.measurement_association.measurement_join_cols": "Measurement-table keys paired positionally with metadata_join_cols.",
    "optional.measurement_association.ordination_methods": "Comma-separated constrained ordination methods to run: cca, rda, and/or dbrda.",
    "optional.grouping_diagnostics.soft_labeling.apply_downstream": "Apply validated proposed labels to downstream metadata and long ASV annotations; false keeps diagnostic proposals separate from analysis labels.",
    "optional.grouping_diagnostics.soft_labeling.exclude_labels": "Labels excluded from soft-label training and predictions; observed excluded labels are retained rather than overwritten.",
    "optional.grouping_diagnostics.soft_labeling.min_cv_balanced_accuracy": "Minimum cross-validated balanced accuracy required before proposed labels can be applied downstream.",
    "optional.asv_mag_link.master_tsv": "Optional precomputed ASV-MAG master-link TSV, used instead of recomputing links from genome/Barrnap inputs.",
    "optional.asv_mag_network.graph_variant": "Network graph to annotate: all or the thresholded subgraph.",
    "optional.asv_mag_network.functional_annotations": "Optional functional-annotation tables to summarize across linked MAGs and network modules.",
    "core.paths.input_dir": "Directory searched for FASTQ files when no manifest is supplied.",
    "core.paths.output_dir": "Public output directory published atomically after a successful run.",
    "core.paths.manifest": "Optional TSV declaring sample IDs and R1/R2 FASTQ paths.",
    "core.paths.runtime_dir": "Persistent private runtime root for work, staging, and environment caches.",
    "core.paths.keep_runtime_dir": "Retain the private runtime after success so resume and targeted reruns remain possible.",
    "core.paths.work_dir": "Optional override for the Nextflow work directory; null derives it from runtime_dir.",
    "core.paths.conda_cache_dir": "Optional override for the per-process Conda cache; null derives it from runtime_dir.",
    "core.resources.threads": "CPU threads requested by each parallel per-sample read-processing task.",
    "core.resources.single_end": "Treat inputs as single-end reads instead of the default paired-end design.",
    "core.filename_patterns.r1_tokens": "Tokens used to recognize read-1 FASTQ filenames during discovery.",
    "core.filename_patterns.r2_tokens": "Read-2 tokens paired positionally with r1_tokens.",
    "core.filename_patterns.ext_patterns": "Regular expressions accepted as FASTQ filename suffixes.",
    "core.filename_patterns.sample_strip_regex": "Regular expression removed from discovered filenames to derive sample IDs.",
    "core.table_filter.min_sample_reads": "General sample-depth cutoff used only when control decontamination is disabled; otherwise min_biological_reads performs biological-only inclusion upstream.",
    "core.table_filter.min_relative_abundance_pct": "Initial per-sample relative-abundance percentage inside FILTER_ASVS; retain an ASV meeting this in any sample. Zero disables abundance exclusion.",
    "core.table_filter.script": "Repository-relative implementation of the initial step inside FILTER_ASVS.",
    "core.filter.max_ee": "Maximum expected errors accepted for a merged read.",
    "core.filter.min_len": "Minimum accepted merged-read length in bases.",
    "core.filter.max_len": "Maximum accepted merged-read length in bases.",
    "core.unoise.min_size": "Minimum dereplicated abundance supplied to UNOISE denoising.",
    "core.swarm.distance": "Sequence-distance setting used for SWARM-style clustering/mapping behavior.",
    "standard.mito.run_mitomaster": "Contact MITOMASTER for mitochondrial evidence; disable for an intentionally offline run.",
    "standard.mito.min_pident": "Minimum BLAST percent identity for local non-target evidence.",
    "standard.mito.min_percov": "Minimum BLAST query coverage percentage for local non-target evidence.",
    "standard.filter_counts.min_relative_abundance_pct": "Per-sample microbial relative-abundance percentage; retain an ASV meeting it in any biological sample after control filtering. General and mock default 0.1%.",
    "standard.filter_counts.min_consensus": "Minimum taxonomy consensus score accepted by final count filtering.",
    "standard.filter_counts.min_group_size": "Minimum number of samples required for a metadata group to participate in group-aware filtering.",
    "standard.filter_counts.exclude_taxa": "Exact rank-qualified taxa removed even when they pass other filters.",
    "optional.voc_correlation.isa_focus_groups": "Sample-type labels defining the focus-group ISA/VOC subset; matching groups may include these labels in singleton or mixed memberships.",
    "optional.voc_correlation.isa_exclude_all_types_from_focus": "Exclude universal sample-type memberships from the focus-group ISA/VOC subset.",
    "optional.voc_correlation.correlation_direction": "Direction retained in the baseline ASV–VOC correlation outputs: positive, negative, or both.",
    "optional.voc_correlation.isa_correlation_direction": "Direction retained in ISA-focused ASV–VOC outputs.",
    "optional.voc_correlation.isa_min_abs_rho": "Minimum absolute Spearman rho required in focused ISA/VOC rows and columns.",
    "optional.voc_correlation.sample_min_abs_z": "Minimum absolute sample VOC z-score used for focused sample heatmaps.",
    "optional.voc_correlation.patient_inference": "Add patient-level relative-abundance permutation correlations and CLR sensitivity analysis; legacy sample correlations remain exploratory.",
    "optional.voc_correlation.patient_permutations": "Seeded permutations for patient correlations and case-status tests; case-status enumeration is exact when all allocations fit this budget. Use 999 for demonstrations, 9999 or more for analysis.",
    "optional.voc_correlation.patient_seed": "Random seed for patient-level VOC permutation tests.",
    "optional.voc_correlation.patient_min_patients": "Minimum patients with cognate VOC measurements to test an ASV–VOC association (at least 3; default 6 is a feasibility gate, not a power guarantee).",
    "optional.voc_correlation.patient_min_nonzero": "Minimum patients with any detected counts for a candidate ASV; insufficient pairs are reported without p-values.",
    "optional.voc_correlation.clr_pseudocount": "Positive count pseudocount added to every supplied ASV before sample-wise centered log-ratios; sensitivity analysis only, default 0.5.",
    "optional.spieceasi.min_relative_abundance_fraction": "Network-input per-sample abundance fraction; 0.001 means 0.1%. Require at least one sample to meet it; force-kept indicators bypass it.",
    "optional.spieceasi.min_prevalence_fraction": "Fraction of network-input samples with a nonzero count; 0.05 means 5%. Force-kept indicators bypass it; zero-variance filtering still applies.",
    "optional.spieceasi.edge_threshold": "Absolute association-weight threshold used to retain reported network edges.",
    "optional.spieceasi.keep_negative": "Retain negative as well as positive inferred network edges.",
    "optional.network_topology.n_null": "Number of seeded degree-preserving null-network draws.",
    "optional.network_topology.skip_null": "Report observed topology without generating the null ensemble.",
    "optional.asv_mag_link.min_pident": "Minimum ASV-to-SSU percent identity accepted as a MAG link.",
    "optional.asv_mag_link.min_qcov": "Minimum ASV query coverage accepted as a MAG link.",
}

COMMON = {
    "enabled": "Whether this module or nested analysis is scheduled.",
    "sample_col": "Column containing sample identifiers.",
    "sample_id_col": "Column containing sample identifiers.",
    "metadata_sample_col": "Sample-identifier column in the metadata table.",
    "patient_col": "Column containing participant/patient identifiers for blocking or pairing.",
    "block_col": "Column defining repeated-measure or permutation blocks.",
    "case_col": "Column containing case/control status.",
    "type_col": "Column containing sample-type labels.",
    "group_col": "Primary grouping column used by this module.",
    "group1_col": "Primary display or analysis grouping column.",
    "group2_col": "Secondary display or analysis grouping column.",
    "color_col": "Metadata column containing or selecting display colors.",
    "count_col": "Column containing ASV counts or transformed count values.",
    "taxon_col": "Column containing taxonomy strings.",
    "feature_col": "Column containing ASV/feature identifiers.",
    "consensus_col": "Column containing taxonomy consensus scores.",
    "formats": "Comma-separated output figure formats.",
    "threads": "CPU threads requested by this module.",
    "ncores": "Parallel worker count requested by this module.",
    "seed": "Random seed used to make stochastic behavior reproducible.",
    "random_state": "Random seed used to make stochastic behavior reproducible.",
    "permutations": "Number of permutations used by the statistical test.",
    "perms": "Number of permutations used by the statistical test.",
    "dpi": "Raster figure resolution in dots per inch.",
    "verbose": "Emit additional module diagnostics.",
    "metadata": "Path to the metadata TSV consumed by this module.",
    "output": "Output filename written by this module.",
    "reference": "Local reference artifact; null permits configured download behavior.",
    "relabel": "Rewrite concatenated FASTA headers so each sequence retains its sample identity.",
    "label_sep": "Separator placed between sample labels and original sequence identifiers.",
    "trim_front_r1": "Number of bases removed from the 5′ end of read 1.",
    "trim_tail_r1": "Number of bases removed from the 3′ end of read 1.",
    "trim_front_r2": "Number of bases removed from the 5′ end of read 2.",
    "trim_tail_r2": "Number of bases removed from the 3′ end of read 2.",
    "trunc_quality": "Quality score used by the read-merging truncation rule.",
    "allow_stagger": "Permit staggered paired-read alignments during merging.",
    "normalize": "Normalization method applied before this analysis or visualization.",
    "transform": "Transformation applied to input abundance values.",
    "transpose": "Transpose the configured matrix orientation before analysis.",
    "figsize": "Figure width and height specification.",
    "style": "Plotting style applied to generated figures.",
    "title": "Title printed on generated figures.",
}


def leaves(value: Any, prefix: str = "") -> Iterator[tuple[str, Any]]:
    """Yield every explicit configurable leaf; lists are single parameters."""
    if isinstance(value, dict):
        if not value:
            yield prefix, value
            return
        for key, child in value.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            yield from leaves(child, path)
    else:
        yield prefix, value


def humanize(key: str) -> str:
    return key.replace("_", " ").replace("asv", "ASV").replace("voc", "VOC")


def describe(path: str, value: Any) -> str:
    if path in EXACT:
        return EXACT[path]
    key = path.rsplit(".", 1)[-1]
    module = path.split(".")[1] if path.startswith(("core.", "standard.", "optional.")) else "workflow"
    label = humanize(key)
    if path.startswith("environments."):
        return f"Conda environment YAML used by the {humanize(key)} process family."
    if key in COMMON:
        return COMMON[key]
    if key.endswith("_col"):
        return f"Input-table column containing {humanize(key[:-4])}."
    if key.endswith("_cols"):
        return f"Input-table columns used for {humanize(key[:-5])}."
    if key.endswith("_palette") or key.endswith("_palettes"):
        return f"Explicit label-to-color mapping used for {label.rsplit(' ', 1)[0]} displays."
    if key.endswith("_order") or key.endswith("_orders"):
        return f"Explicit category order used for {label.rsplit(' ', 1)[0]} analysis and display."
    if key.endswith("_dir") or key.endswith("_subdir") or key == "sub_dir":
        return f"Directory or published subdirectory used for {label}."
    if key.endswith("_file") or key.endswith("_table") or key.endswith("_input"):
        return f"Path to, or configured name of, the {label} input."
    if key.endswith("_fasta"):
        return f"Path to the {humanize(key[:-6])} FASTA reference; null means it is not supplied."
    if key.endswith("_db"):
        return f"BLAST database prefix used for {humanize(key[:-3])} screening."
    if key.endswith("_url"):
        return f"Download URL used when the local {humanize(key[:-4])} is unavailable."
    if key.endswith("_filename"):
        return f"Local filename assigned to the retrieved {humanize(key[:-9])} artifact."
    if key.endswith("_output") or key.endswith("_output_dir") or key.endswith("_prefix"):
        return f"Output name, prefix, or destination used for {label}."
    if key.startswith("min_"):
        return f"Minimum accepted {humanize(key[4:])} for the {humanize(module)} module."
    if key.startswith("max_"):
        return f"Maximum accepted {humanize(key[4:])} for the {humanize(module)} module."
    if key.startswith("skip_"):
        return f"Skip {humanize(key[5:])} when true."
    if key.startswith("run_"):
        return f"Run {humanize(key[4:])} when true."
    if key.startswith("keep_"):
        return f"Retain {humanize(key[5:])} when true."
    if key.startswith("exclude_"):
        return f"Values or groups excluded according to {humanize(key[8:])}."
    if key.startswith("force_"):
        return f"Force {humanize(key[6:])} instead of accepting a reusable cached product."
    if key.endswith("_threshold") or key.endswith("_alpha"):
        return f"Statistical or filtering cutoff for {label} in the {humanize(module)} module."
    if key.endswith("_size") or key.endswith("_n") or key.startswith("n_"):
        return f"Configured number or size for {label} in the {humanize(module)} module."
    if isinstance(value, bool):
        return f"Enable or disable {label} behavior in the {humanize(module)} module."
    if isinstance(value, list):
        return f"Ordered values used for {label} by the {humanize(module)} module."
    if isinstance(value, (int, float)):
        return f"Numeric setting for {label} in the {humanize(module)} module."
    return f"Value selecting or naming {label} for the {humanize(module)} module."


def render(value: Any) -> str:
    if value is None:
        shown = "null"
    elif isinstance(value, bool):
        shown = str(value).lower()
    elif isinstance(value, (list, dict)):
        shown = json.dumps(value, separators=(",", ":"))
    else:
        shown = str(value)
    return shown.replace("|", "\\|").replace("\n", " ")


def main() -> None:
    config = yaml.safe_load(TEMPLATE.read_text())
    lines = [
        "# Complete ASPIRE Configuration Parameter Catalogue",
        "",
        "This file documents every explicit parameter in the canonical",
        "{download}`asv_pipeline_nextflow.yml <../asv_pipeline_nextflow.yml>` template. It is generated",
        "from that template; lists are documented as one parameter and their complete template",
        "value is shown. Paths are resolved relative to the YAML file unless absolute.",
        "",
        "The template value is a starting point, not a universal scientific recommendation.",
        "Study-specific thresholds, metadata columns, references, and optional modules must be",
        "chosen and reported for the study design. See [Configuration Reference](CONFIGURATION.md)",
        "for dependencies and interpretation and [Process Reference](PROCESS_REFERENCE.md) for I/O.",
        "",
    ]
    for tier, tier_config in config.items():
        lines.extend([f"## `{tier}`", ""])
        if tier == "environments":
            sections = {"process environments": tier_config}
        else:
            sections = tier_config
        for section, section_config in sections.items():
            lines.extend([
                f"### `{tier}.{section}`" if tier != "environments" else "### Process environment definitions",
                "",
            ])
            if tier == "optional" and section == "measurement_association":
                lines.extend(["This general measurement-association extension is outside the publication workflow. Use `optional.voc_correlation` for ASV associations with Volatile Organic Compounds (VOCs).", ""])
            lines.extend([
                "| Parameter | Type | Template value | Definition |",
                "|---|---|---|---|",
            ])
            prefix = f"{tier}.{section}" if tier != "environments" else tier
            for path, value in leaves(section_config, prefix):
                value_type = "null" if value is None else "list" if isinstance(value, list) else "mapping" if isinstance(value, dict) else type(value).__name__
                lines.append(f"| `{path}` | {value_type} | `{render(value)}` | {describe(path, value)} |")
            lines.append("")
    OUTPUT.write_text("\n".join(lines).rstrip() + "\n")


if __name__ == "__main__":
    main()
