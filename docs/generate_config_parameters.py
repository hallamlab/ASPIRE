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
    "core.table_filter.min_sample_sum": "Minimum total reads required for a sample to survive the early technical count filter.",
    "core.table_filter.min_asv_sum": "Minimum study-wide count required for an ASV to survive the early technical count filter.",
    "core.table_filter.script": "Repository-relative implementation used for early count-table filtering.",
    "core.filter.max_ee": "Maximum expected errors accepted for a merged read.",
    "core.filter.min_len": "Minimum accepted merged-read length in bases.",
    "core.filter.max_len": "Maximum accepted merged-read length in bases.",
    "core.unoise.min_size": "Minimum dereplicated abundance supplied to UNOISE denoising.",
    "core.swarm.distance": "Sequence-distance setting used for SWARM-style clustering/mapping behavior.",
    "standard.mito.run_mitomaster": "Contact MITOMASTER for mitochondrial evidence; disable for an intentionally offline run.",
    "standard.mito.min_pident": "Minimum BLAST percent identity for local non-target evidence.",
    "standard.mito.min_percov": "Minimum BLAST query coverage percentage for local non-target evidence.",
    "standard.filter_counts.abundance_threshold": "Minimum relative-abundance percentage used by the final abundance filter; 0.5 means 0.5%.",
    "standard.filter_counts.min_consensus": "Minimum taxonomy consensus score accepted by final count filtering.",
    "standard.filter_counts.min_group_size": "Minimum number of samples required for a metadata group to participate in group-aware filtering.",
    "standard.filter_counts.exclude_taxa": "Exact rank-qualified taxa removed even when they pass other filters.",
    "optional.three_tier_decontam.pooled_threshold": "Decontam score cutoff for the pooled control-based tier.",
    "optional.three_tier_decontam.within_type_threshold": "Decontam score cutoff applied within each configured sample type.",
    "optional.three_tier_decontam.aggressive_threshold": "Secondary, more permissive contaminant score cutoff used by the configured combination rule.",
    "optional.three_tier_decontam.combine_mode": "Rule used to combine pooled and within-type contaminant evidence.",
    "optional.three_tier_decontam.biological_plausibility": "Apply the taxonomy/biological-plausibility tier before producing downstream tables.",
    "optional.voc_correlation.correlation_direction": "Direction retained in the baseline ASV–VOC correlation outputs: positive, negative, or both.",
    "optional.voc_correlation.isa_correlation_direction": "Direction retained in ISA-focused ASV–VOC outputs.",
    "optional.voc_correlation.isa_min_abs_rho": "Minimum absolute Spearman rho required in focused ISA/VOC rows and columns.",
    "optional.voc_correlation.sample_min_abs_z": "Minimum absolute sample VOC z-score used for focused sample heatmaps.",
    "optional.voc_correlation.patient_inference": "Add patient-level relative-abundance permutation correlations and CLR sensitivity analysis; legacy sample correlations remain exploratory.",
    "optional.voc_correlation.patient_permutations": "Seeded permutations for patient correlations and case-status tests; case-status enumeration is exact when all allocations fit this budget. Use 999 for demonstrations, 9999 or more for analysis.",
    "optional.voc_correlation.patient_seed": "Random seed for patient-level VOC permutation tests.",
    "optional.voc_correlation.patient_min_patients": "Minimum patients with cognate measurements to test an ASV–VOC association (at least 3; default 6 is a feasibility gate, not a power guarantee).",
    "optional.voc_correlation.patient_min_nonzero": "Minimum patients with any detected counts for a candidate ASV; insufficient pairs are reported without p-values.",
    "optional.voc_correlation.clr_pseudocount": "Positive count pseudocount added to every supplied ASV before sample-wise centered log-ratios; sensitivity analysis only, default 0.5.",
    "optional.spieceasi.min_rel_abund": "Minimum relative abundance required for an ASV to enter network inference.",
    "optional.spieceasi.min_prevalence": "Minimum prevalence required for an ASV to enter network inference.",
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
        "[`asv_pipeline_nextflow.yml`](../asv_pipeline_nextflow.yml) template. It is generated",
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
