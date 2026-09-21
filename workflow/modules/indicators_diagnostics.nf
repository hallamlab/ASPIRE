process INDICSPECIES {
    cpus pipelineThreads
    conda "${indicspeciesCondaEnvPath}"

    input:
    path(metadata_table)
    path(asv_counts)

    output:
    path("indicspecies_group1_summary.tsv"), emit: group1_summary
    path("indicspecies_group2_summary.tsv"), emit: group2_summary
    path("indicspecies_group1_results.tsv"), emit: group1_results
    path("indicspecies_group2_results.tsv"), emit: group2_results
    path("indicspecies_tables/*.tsv"), emit: all_tables
    path("indicspecies.done"), emit: done

    when:
    indicspeciesEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def indicspeciesGroupColsArg = indicspeciesGroupCols.join(',')
    def indicspeciesBlockArg = indicspeciesBlockCol ? """  --block-col "${indicspeciesBlockCol}" \\\n""" : ''
    def indicspeciesBlockedArg = indicspeciesBlockedCols ? """  --blocked-cols "${indicspeciesBlockedCols}" \\\n""" : ''
    def indicspeciesStratifiedArg = indicspeciesStratifiedSpecsArg ? """  --stratified-isa "${indicspeciesStratifiedSpecsArg}" \\\n""" : ''
    def isaSummarySuffix = '_indicator_species_summary.tsv'
    def isaResultsSuffix = '_indicator_species_results.tsv'
    def group1SummaryPath = "${indicspeciesOutputDirAbs}/${indicspeciesGroup1}${isaSummarySuffix}"
    def group2SummaryPath = "${indicspeciesOutputDirAbs}/${indicspeciesGroup2}${isaSummarySuffix}"
    def group1ResultsPath = "${indicspeciesOutputDirAbs}/${indicspeciesGroup1}${isaResultsSuffix}"
    def group2ResultsPath = "${indicspeciesOutputDirAbs}/${indicspeciesGroup2}${isaResultsSuffix}"
    """
set -euo pipefail
ISA_STAGE_ROOT="\$PWD/isa_stage_root"
mkdir -p "\${ISA_STAGE_ROOT}"

Rscript "${indicspeciesScriptPath}" \\
  --asv "${asv_counts}" \\
  --meta "${metadata_table}" \\
  --sample-col "${indicspeciesSampleCol}" \\
  --group-cols "${indicspeciesGroupColsArg}" \\
${indicspeciesBlockArg}${indicspeciesBlockedArg}${indicspeciesStratifiedArg}  --perms ${indicspeciesPerms} \\
  --seed ${indicspeciesSeed} \\
  --q-threshold ${indicspeciesQThreshold} \\
  --min-n ${indicspeciesMinN} \\
  --outdir "\${ISA_STAGE_ROOT}"

"\${CONDA_PREFIX}/bin/python" - <<'PY'
from pathlib import Path
import shutil
import sys

staged = Path("isa_stage_root/indicspecies")
final = Path("${indicspeciesOutputDirAbs}")
if not staged.is_dir():
    print(f"Missing staged indicspecies directory: {staged}", file=sys.stderr)
    raise SystemExit(1)
required_names = [
    Path("${group1SummaryPath}").name,
    Path("${group2SummaryPath}").name,
    Path("${group1ResultsPath}").name,
    Path("${group2ResultsPath}").name,
]
missing_staged = [name for name in required_names if not (staged / name).is_file()]
if missing_staged:
    for name in missing_staged:
        print(f"Missing expected staged indicspecies output: {staged / name}", file=sys.stderr)
    raise SystemExit(1)
if final.exists():
    shutil.rmtree(final)
final.parent.mkdir(parents=True, exist_ok=True)
shutil.move(str(staged), str(final))

group1_summary = Path("${group1SummaryPath}")
group2_summary = Path("${group2SummaryPath}")
group1_results = Path("${group1ResultsPath}")
group2_results = Path("${group2ResultsPath}")
required = [group1_summary, group2_summary, group1_results, group2_results]

missing = [str(p) for p in required if not p.is_file()]
if missing:
    for p in missing:
        print(f"Missing expected indicspecies output: {p}", file=sys.stderr)
    raise SystemExit(1)

tables_dir = Path("indicspecies_tables")
tables_dir.mkdir(parents=True, exist_ok=True)
all_tables = sorted(Path("${indicspeciesOutputDirAbs}").glob("*_indicator_species*.tsv"))
if not all_tables:
    print("No indicspecies tables were generated in ${indicspeciesOutputDirAbs}", file=sys.stderr)
    raise SystemExit(1)

for src in all_tables:
    dst = tables_dir / src.name
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src.resolve())

link_map = {
    "indicspecies_group1_summary.tsv": group1_summary,
    "indicspecies_group2_summary.tsv": group2_summary,
    "indicspecies_group1_results.tsv": group1_results,
    "indicspecies_group2_results.tsv": group2_results,
}
for dst_name, src in link_map.items():
    dst = Path(dst_name)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src.resolve())
PY

touch indicspecies.done
"""
}
process INDICSPECIES_PLOTS {
    cpus pipelineThreads
    conda "${indicspeciesCondaEnvPath}"

    input:
    path(metadata_table)
    path(indicspecies_tables)

    output:
    path("indicspecies_plots.done"), emit: done

    when:
    indicspeciesEnabled && indicspeciesPlotEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def plotVennPath = indicspeciesPlotVennPath ?: ''
    def plotTaxPath = indicspeciesPlotTaxonomyPath ?: ''
    def plotPairsMode = indicspeciesPlotPairsMode?.toString()?.trim()?.toLowerCase() ?: 'all'
    """
set -euo pipefail
mkdir -p "${indicspeciesPlotOutputDirAbs}"

"\${CONDA_PREFIX}/bin/python" - <<'PY'
from pathlib import Path
import itertools
import json
import re
import shutil
import subprocess
import sys

plot_script = Path("${plotIndicspeciesScriptPath}")
final_out_root = Path("${indicspeciesPlotOutputDirAbs}")
out_root = Path("indicspecies_plots_staged")
out_root.mkdir(parents=True, exist_ok=True)
pairs_mode = "${plotPairsMode}"
plot_tax = Path("${plotTaxPath}") if "${plotTaxPath}" else None
plot_venn = Path("${plotVennPath}") if "${plotVennPath}" else None
metadata_path = Path("${metadata_table}") if "${metadata_table}" else None
preferred_group1 = "${indicspeciesGroup1}".strip()
preferred_group2 = "${indicspeciesGroup2}".strip()
metadata_color_col = "${indicspeciesColorCol}".strip()
palette_cfg = json.loads(r'''${indicspeciesGroupPaletteJson}''')
order_cfg = json.loads(r'''${indicspeciesGroupOrderJson}''')
focus_cfg = json.loads(r'''${indicspeciesFocusLabelJson}''')
label_focused_asvs = ${indicspeciesLabelFocusedAsvs ? 'True' : 'False'}
configured_groups = json.loads(r'''${groovy.json.JsonOutput.toJson(indicspeciesGroupCols)}''')
selected_suffix = "_indicator_species_summary.tsv"

summary_files = [Path(f"{group}{selected_suffix}") for group in configured_groups]
missing_summaries = [path.name for path in summary_files if not path.is_file()]
if missing_summaries:
    print("[STOP] Missing configured ISA summaries: " + ", ".join(missing_summaries), file=sys.stderr)
    raise SystemExit(1)
if len(summary_files) < 2:
    print(
        f"[w] Need at least 2 indicspecies summary tables to build ISA plots; found {len(summary_files)}. Skipping.",
        file=sys.stderr,
    )
    raise SystemExit(0)

def clean_group_name(raw_name: str) -> str:
    name = raw_name
    for ending in ("_summary.tsv", "_results.tsv"):
        if name.endswith(ending):
            name = name[: -len(ending)]
            break
    for suffix in ("_indicator_species_DULEG", "_indicator_species"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name

def orient_pair(a: Path, b: Path) -> tuple[Path, Path]:
    an = clean_group_name(a.name)
    bn = clean_group_name(b.name)
    if an == preferred_group1 and bn == preferred_group2:
        return a, b
    if an == preferred_group2 and bn == preferred_group1:
        return b, a
    if an == preferred_group1:
        return a, b
    if bn == preferred_group1:
        return b, a
    return a, b

if pairs_mode == "first_vs_rest":
    raw_pairs = [(summary_files[0], f) for f in summary_files[1:]]
else:
    raw_pairs = list(itertools.combinations(summary_files, 2))
pair_iter = [orient_pair(a, b) for a, b in raw_pairs]

# Record the bounded plot plan before rendering. Stratified ISA tables are
# intentionally excluded: crossing every stratum with every unrelated ISA is
# neither interpretable nor computationally bounded.
(out_root / "isa_plot_plan.tsv").write_text(
    chr(9).join(("group1", "group2", "group1_file", "group2_file")) + chr(10) +
    "".join(
        chr(9).join((clean_group_name(a.name), clean_group_name(b.name), a.name, b.name)) + chr(10)
        for a, b in pair_iter
    ),
    encoding="utf-8",
)

def has_col(path: Path, col: str) -> bool:
    if not col:
        return False
    try:
        header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    except Exception:
        return False
    return col in header

def has_meta_col(path: Path | None, col: str) -> bool:
    if path is None or not path.is_file() or not col:
        return False
    return has_col(path, col)

def order_string(group_name: str) -> str:
    raw = order_cfg.get(group_name, [])
    if isinstance(raw, list):
        return ",".join(str(x).strip() for x in raw if str(x).strip())
    if raw:
        return str(raw).strip()
    return ""

for g1_file, g2_file in pair_iter:
    g1_name = clean_group_name(g1_file.name)
    g2_name = clean_group_name(g2_file.name)
    if g1_name == g2_name:
        print(
            f"[w] Skipping same-group ISA pair: {g1_file.name} vs {g2_file.name}",
            file=sys.stderr,
        )
        continue
    pair_slug = re.sub(r"[^0-9A-Za-z._-]", "_", f"{g1_name}__{g2_name}")
    pair_out = out_root / pair_slug
    pair_out.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(plot_script),
        "--group1-results",
        str(g1_file),
        "--group2-results",
        str(g2_file),
        "--group1-name",
        g1_name,
        "--group2-name",
        g2_name,
        "--outdir",
        str(pair_out),
    ]
    if metadata_path and metadata_path.is_file():
        cmd.extend(["--metadata", str(metadata_path)])
        if has_meta_col(metadata_path, g1_name):
            cmd.extend(["--group1-meta-label-col", g1_name])
        if has_meta_col(metadata_path, g2_name):
            cmd.extend(["--group2-meta-label-col", g2_name])
        if g1_name == preferred_group1 and metadata_color_col and has_meta_col(metadata_path, metadata_color_col):
            cmd.extend(["--group1-meta-color-col", metadata_color_col])

    if has_col(g1_file, g1_name):
        cmd.extend(["--group1-label-col", g1_name])
    if has_col(g1_file, f"{g1_name}_Color"):
        cmd.extend(["--group1-color-col", f"{g1_name}_Color"])
    elif has_col(g1_file, "Color"):
        cmd.extend(["--group1-color-col", "Color"])

    if has_col(g2_file, g2_name):
        cmd.extend(["--group2-label-col", g2_name])
    if has_col(g2_file, f"{g2_name}_Color"):
        cmd.extend(["--group2-color-col", f"{g2_name}_Color"])
    elif has_col(g2_file, "Color"):
        cmd.extend(["--group2-color-col", "Color"])

    if has_col(g2_file, f"{g2_name}_Marker"):
        cmd.extend(["--group2-marker-col", f"{g2_name}_Marker"])
    if palette_cfg.get(g1_name):
        cmd.extend(["--group1-palette", str(palette_cfg[g1_name])])
    if palette_cfg.get(g2_name):
        cmd.extend(["--group2-palette", str(palette_cfg[g2_name])])
    group1_order_cfg = order_string(g1_name)
    group2_order_cfg = order_string(g2_name)
    if group1_order_cfg:
        cmd.extend(["--group1-order", group1_order_cfg])
    if group2_order_cfg:
        cmd.extend(["--group2-order", group2_order_cfg])
    if focus_cfg.get(g1_name):
        cmd.extend(["--focus-group1-label", str(focus_cfg[g1_name])])
    if focus_cfg.get(g2_name):
        cmd.extend(["--focus-group2-label", str(focus_cfg[g2_name])])
    if label_focused_asvs:
        cmd.append("--label-focused-asvs")
    if plot_tax and plot_tax.is_file():
        cmd.extend(["--taxonomy", str(plot_tax)])
    if plot_venn and plot_venn.is_file():
        cmd.extend(["--venn", str(plot_venn)])

    subprocess.run(cmd, check=True)

# Replace the prior plot tree only after every requested figure succeeds. This
# removes obsolete combinatorial plot directories from earlier implementations
# without risking partial publication output when a renderer fails.
if final_out_root.exists():
    shutil.rmtree(final_out_root)
final_out_root.parent.mkdir(parents=True, exist_ok=True)
shutil.move(str(out_root), str(final_out_root))
PY

touch indicspecies_plots.done
"""
}
process INDICSPECIES_ALIGNED_PLOTS {
    cpus pipelineThreads
    conda "${indicspeciesCondaEnvPath}"

    input:
    path(indicspecies_tables)

    output:
    path("indicspecies_aligned.done"), emit: done

    when:
    indicspeciesEnabled && indicspeciesAlignedEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
mkdir -p "${indicspeciesAlignedOutputDirAbs}"
mkdir -p aligned_indicspecies_input

for src in *.tsv; do
  [[ -f "\${src}" ]] || continue
  ln -sf "\$(realpath "\${src}")" "aligned_indicspecies_input/\$(basename "\${src}")"
done

"\${CONDA_PREFIX}/bin/python" "${plotIndicspeciesAlignedScriptPath}" \\
  --indicspecies-dir aligned_indicspecies_input \\
  --outdir "${indicspeciesAlignedOutputDirAbs}" \\
  --alpha ${indicspeciesAlignedAlpha} \\
  --min-stat ${indicspeciesAlignedMinStat} \\
  --top-n ${indicspeciesAlignedTopN}

touch indicspecies_aligned.done
"""
}
process VOC_CORRELATION {
    cpus pipelineThreads
    conda "${vocCorrelationCondaEnvPath}"

    input:
    path(asv_meta_table)
    path(asv_counts)
    path(indicspecies_tables)

    output:
    path("voc_correlation.done"), emit: done

    when:
    vocCorrelationEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def vocColsArgs = vocCorrelationVocCols.collect { col -> """  --voc-col "${col}" \\\n""" }.join('')
    def legacySubsetArg = vocCorrelationUseLegacySubset ? "  --use-legacy-voc-subset \\\n" : ''
"""
set -euo pipefail
mkdir -p "${vocCorrelationOutputDirAbs}"
echo "plot_voc_corr.py md5: ${plotVocCorrScriptHash}"

python3 "${plotVocCorrScriptPath}" \\
  --asv-meta "${asv_meta_table}" \\
  --asv-counts "${asv_counts}" \\
  --voc "${vocCorrelationVocTablePath}" \\
  --outdir "${vocCorrelationOutputDirAbs}" \\
  --metadata-sample-col "${vocCorrelationMetadataSampleCol}" \\
  --type-col "${vocCorrelationTypeCol}" \\
  --patient-col "${vocCorrelationPatientCol}" \\
  --case-col "${vocCorrelationCaseCol}" \\
  --sample-types "${vocCorrelationSampleTypes}" \\
  --voc-sample-col "${vocCorrelationVocSampleCol}" \\
  --sample-id-mode "${vocCorrelationSampleIdMode}" \\
${legacySubsetArg}${vocColsArgs}  --spieceasi-min-rel-abund ${spieceasiMinRelAbund} \\
  --spieceasi-min-prevalence ${spieceasiMinPrevalence} \\
  --spieceasi-remove-zero-var ${spieceasiRemoveZeroVar} \\
  --correlation-direction "${vocCorrelationDirection}" \\
  --case-palette "${vocCorrelationCasePalette}" \\
  --isa-palette "${vocCorrelationIsaPalette}" \\
  --isa-q-threshold ${indicspeciesQThreshold} \\
  --indicspecies-glob "*_indicator_species*.tsv"

touch voc_correlation.done
"""
}
process MEASUREMENT_ASSOCIATION {
    cpus pipelineThreads
    conda "${measurementAssociationCondaEnvPath}"

    input:
    path(asv_meta_table)
    path(metadata_table)
    path(asv_counts)
    path(asv_subset_audit)

    output:
    path("measurement_association.done"), emit: done

    when:
    measurementAssociationEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def measurementTableArg = measurementAssociationTablePath ? """  --measurement-table "${measurementAssociationTablePath}" \\\n""" : ''
    def measurementColsArg = measurementAssociationCols ? """  --measurement-cols "${measurementAssociationCols.join('|')}" \\\n""" : ''
    def excludeColsArg = measurementAssociationExcludeCols ? """  --exclude-cols "${measurementAssociationExcludeCols.join('|')}" \\\n""" : ''
    def metadataJoinArg = measurementAssociationMetadataJoinCols ? measurementAssociationMetadataJoinCols.join(',') : ''
    def measurementJoinArg = measurementAssociationMeasurementJoinCols ? measurementAssociationMeasurementJoinCols.join(',') : ''
    def subsetAuditArg = measurementAssociationSubsetSource == 'spieceasi_standard_filter' ? """  --asv-subset-audit "${asv_subset_audit}" \\\n""" : ''
"""
set -euo pipefail
mkdir -p "${measurementAssociationOutputDirAbs}"
echo "measurement_association.py md5: ${measurementAssociationScriptHash}"
echo "run_measurement_association.R md5: ${measurementAssociationRScriptHash}"

"\${CONDA_PREFIX}/bin/python" "${measurementAssociationScriptPath}" \\
  --asv-meta "${asv_meta_table}" \\
  --metadata "${metadata_table}" \\
  --asv-counts "${asv_counts}" \\
${measurementTableArg}  --outdir "${measurementAssociationOutputDirAbs}" \\
  --r-script "${measurementAssociationRScriptPath}" \\
  --sample-col "${measurementAssociationSampleCol}" \\
  --asv-id-col "${measurementAssociationAsvIdCol}" \\
  --measurement-sample-col "${measurementAssociationMeasurementSampleCol}" \\
  --metadata-join-cols "${metadataJoinArg}" \\
  --measurement-join-cols "${measurementJoinArg}" \\
  --measurement-aliases-json '${measurementAssociationAliasesJson}' \\
${measurementColsArg}${excludeColsArg}  --group-col "${measurementAssociationGroupCol}" \\
  --group-palette "${measurementAssociationGroupPalette}" \\
  --asv-subset-source "${measurementAssociationSubsetSource}" \\
${subsetAuditArg}\
  --max-asvs ${measurementAssociationMaxAsvs} \\
  --min-total ${measurementAssociationMinTotal} \\
  --min-prevalence ${measurementAssociationMinPrevalence} \\
  --top-correlations ${measurementAssociationTopCorrelations} \\
  --correlation-direction "${measurementAssociationDirection}" \\
  --ordination-methods "${measurementAssociationMethods}" \\
  --ordination-measurement-cols "${measurementAssociationOrdinationCols}" \\
  --permutations ${measurementAssociationPermutations} \\
  --top-vectors ${measurementAssociationTopVectors} \\
  --formats "${measurementAssociationFormats}"

touch measurement_association.done
"""
}

process MODULE_MEASUREMENT_ASSOCIATION {
    cpus 1
    conda "${measurementAssociationCondaEnvPath}"

    input:
    path(asv_counts)
    path(prevalence_counts)
    path(metadata_table)
    path(modules_all)
    path(network_graph)
    path(spieceasi_layout)
    path(node_features)
    path(taxonomy_table)
    path(measurement_done)
    val(network_metric_top_n)

    output:
    path("module_measurement_association.done"), emit: done

    when:
    moduleMeasurementAssociationEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def pcaArgs = moduleMeasurementPcaScores ? """  --pca-scores "${moduleMeasurementPcaScores}" \\
  --pca-loadings "${moduleMeasurementPcaLoadings}" \\
  --pca-explained "${moduleMeasurementPcaExplained}" \\
  --hybrid-assignments "${moduleMeasurementHybridAssignments}" \\
  --hybrid-centroids "${moduleMeasurementHybridCentroids}" \\
""" : ''
    """
set -euo pipefail
echo "module_measurement_association.py md5: ${moduleMeasurementAssociationScriptHash}"

"\${CONDA_PREFIX}/bin/python" "${moduleMeasurementAssociationScriptPath}" \\
  --asv-counts "${asv_counts}" \\
  --prevalence-counts "${prevalence_counts}" \\
  --metadata "${metadata_table}" \\
  --modules "${modules_all}" \\
  --network-graph "${network_graph}" \\
  --spieceasi-layout "${spieceasi_layout}" \\
  --node-features "${node_features}" \\
  --taxonomy "${taxonomy_table}" \\
  --anchor-top-n ${asvMagNetworkAnchorTopN} \\
  --network-metric-top-n ${network_metric_top_n} \\
  --measurement-matrix "${measurementAssociationOutputDirAbs}/tables/measurement_matrix.samples_by_measurement.tsv" \\
  --asv-correlations "${measurementAssociationOutputDirAbs}/tables/asv_measurement_spearman_long.tsv" \\
  --measurement-selection-audit "${measurementAssociationOutputDirAbs}/tables/measurement_selection_audit.tsv" \\
  --outdir "${measurementAssociationOutputDirAbs}" \\
  --sample-col "${measurementAssociationSampleCol}" \\
  --cruise-col "${moduleMeasurementCruiseCol}" \\
  --depth-col "${moduleMeasurementDepthCol}" \\
  --sparse-measurements "${moduleMeasurementSparseCols}" \\
  --q-threshold ${moduleMeasurementQThreshold} \\
${pcaArgs}  --formats "${measurementAssociationFormats}"

touch module_measurement_association.done
"""
}

process TITAN_PREPARE {
    cpus 1
    conda "${titanCondaEnvPath}"

    input:
    path(asv_counts)
    path(spieceasi_filter_audit)
    path(taxonomy_table)
    path(measurement_done)

    output:
    path("titan_prepared"), emit: prepared
    path("titan_prepare.done"), emit: done

    when:
    titanEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def variablesArg = titanVariables ? """  --variables "${titanVariables}" \\\n""" : ''
    def transposeArg = titanTranspose ? '--transpose' : '--no-transpose'
    """
set -euo pipefail
echo "prepare_titan.py md5: ${titanPrepareScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"

"\${CONDA_PREFIX}/bin/python" "${titanPrepareScriptPath}" \\
  --counts "${asv_counts}" \\
  --filter-audit "${spieceasi_filter_audit}" \\
  --measurement-matrix "${measurementAssociationOutputDirAbs}/tables/measurement_matrix.samples_by_measurement.tsv" \\
  --taxonomy "${taxonomy_table}" \\
  --outdir titan_prepared \\
  --sample-col "${measurementAssociationSampleCol}" \\
  --asv-id-col "${measurementAssociationAsvIdCol}" \\
${variablesArg}  ${transposeArg} \\
  --min-split ${titanMinSplit} \\
  --minimum-samples ${titanMinimumSamples} \\
  --minimum-occurrence ${titanMinimumOccurrence} \\
  --minimum-prevalence ${titanMinimumPrevalence} \\
  --minimum-mean-relative-abundance ${titanMinimumMeanRelativeAbundance}

touch titan_prepare.done
"""
}

process TITAN_INSTALL {
    cpus 1
    conda "${titanCondaEnvPath}"

    output:
    path("titan_r_library"), emit: library
    path("titan2_installation.tsv"), emit: manifest

    when:
    titanEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
echo "install_titan2.R md5: ${titanInstallScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"

Rscript "${titanInstallScriptPath}" \\
  --library titan_r_library \\
  --repository "${titanSourceRepository}" \\
  --revision "${titanSourceRevision}" \\
  --version "${titanExpectedVersion}" \\
  --manifest titan2_installation.tsv
"""
}

process TITAN_ANALYSIS {
    cpus titanNcpus
    conda "${titanCondaEnvPath}"

    input:
    path(prepared)
    path(titan_library)

    output:
    path("titan_raw"), emit: raw
    path("titan_analysis.done"), emit: done

    when:
    titanEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def imaxArg = titanImax ? 'TRUE' : 'FALSE'
    def ivTotalArg = titanIvTotal ? 'TRUE' : 'FALSE'
    def memoryArg = titanMemory ? 'TRUE' : 'FALSE'
    """
set -euo pipefail
echo "run_titan.R md5: ${titanAnalysisScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"

# TITAN2 creates one R worker per configured CPU. Keep native numerical
# libraries single-threaded inside each worker to prevent nested oversubscription.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

R_LIBS_USER="${titan_library}" Rscript "${titanAnalysisScriptPath}" \\
  --prepared-dir "${prepared}" \\
  --outdir titan_raw \\
  --min-split ${titanMinSplit} \\
  --permutations ${titanPermutations} \\
  --bootstrap-count ${titanBootstrapCount} \\
  --seed ${titanSeed} \\
  --imax ${imaxArg} \\
  --iv-total ${ivTotalArg} \\
  --purity-cutoff ${titanPurityCutoff} \\
  --reliability-cutoff ${titanReliabilityCutoff} \\
  --ncpus ${titanNcpus} \\
  --memory ${memoryArg}

touch titan_analysis.done
"""
}

process TITAN_COLLECT {
    cpus 1
    conda "${titanCondaEnvPath}"

    input:
    path(prepared)
    path(raw)
    path(installation_manifest)

    output:
    path("titan_collect.done"), emit: done
    path("titan_taxon_results.tsv"), emit: taxon_results

    when:
    titanEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
echo "collect_titan.py md5: ${titanCollectScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"
rm -rf titan_staged
"\${CONDA_PREFIX}/bin/python" "${titanCollectScriptPath}" \\
  --prepared-dir "${prepared}" \\
  --raw-dir "${raw}" \\
  --installation-manifest "${installation_manifest}" \\
  --outdir titan_staged

"\${CONDA_PREFIX}/bin/python" -c 'from pathlib import Path; import shutil; src=Path("titan_staged"); dst=Path("${titanOutputDirAbs}"); shutil.rmtree(dst, ignore_errors=True); dst.parent.mkdir(parents=True, exist_ok=True); shutil.copytree(src, dst)'
ln -sf "${titanOutputDirAbs}/tables/titan_taxon_results.tsv" titan_taxon_results.tsv
touch titan_collect.done
"""
}

process TITAN_PLOTS {
    cpus 1
    conda "${titanCondaEnvPath}"

    input:
    path(collect_done)

    output:
    path("titan_plots.done"), emit: done

    when:
    titanEnabled && titanPlotEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
echo "plot_titan.py md5: ${titanPlotScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"
mkdir -p "${titanOutputDirAbs}/plots"
"\${CONDA_PREFIX}/bin/python" "${titanPlotScriptPath}" \\
  --tables-dir "${titanOutputDirAbs}/tables" \\
  --outdir "${titanOutputDirAbs}/plots" \\
  --formats "${titanFormats}" \\
  --ranked-top-n ${titanRankedTopN}
touch titan_plots.done
"""
}

process MICROBIAL_COMPARTMENT_PREPARE {
    cpus 1
    conda "${microbialCompartmentCondaEnvPath}"

    input:
    path(asv_counts)
    path(spieceasi_filter_audit)
    path(metadata_table)

    output:
    path("microbial_compartment_prepared"), emit: prepared
    path("microbial_compartment_prepare.done"), emit: done

    when:
    microbialCompartmentEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def transposeArg = microbialCompartmentTranspose ? '--transpose' : '--no-transpose'
    """
set -euo pipefail
echo "prepare_microbial_compartments.py md5: ${microbialCompartmentPrepareScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"

"\${CONDA_PREFIX}/bin/python" "${microbialCompartmentPrepareScriptPath}" \\
  --counts "${asv_counts}" \\
  --spieceasi-filter-audit "${spieceasi_filter_audit}" \\
  --metadata "${metadata_table}" \\
  --outdir microbial_compartment_prepared \\
  --sample-col "${microbialCompartmentSampleCol}" \\
  --asv-id-col "${measurementAssociationAsvIdCol}" \\
  --stability-block-col "${microbialCompartmentStabilityBlockCol}" \\
  --stability-stratum-col "${microbialCompartmentStabilityStratumCol}" \\
  ${transposeArg} \\
  --minimum-prevalence ${microbialCompartmentMinimumPrevalence} \\
  --minimum-max-relative-abundance ${microbialCompartmentMinimumMaxRelativeAbundance} \\
  --zero-replacement "${microbialCompartmentZeroReplacement}" \\
  --multiplicative-delta ${microbialCompartmentMultiplicativeDelta} \\
  --pseudocount ${microbialCompartmentPseudocount}

touch microbial_compartment_prepare.done
"""
}

process MICROBIAL_COMPARTMENT_INFERENCE {
    cpus microbialCompartmentNcpus
    conda "${microbialCompartmentCondaEnvPath}"

    input:
    path(prepared)

    output:
    path("microbial_compartment_inference"), emit: inference
    path("microbial_compartment_inference.done"), emit: done

    when:
    microbialCompartmentEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def hierarchicalArg = microbialCompartmentHierarchicalEnabled ? 'TRUE' : 'FALSE'
    def predictionStrengthArg = microbialCompartmentPredictionStrengthEnabled ? 'TRUE' : 'FALSE'
    def blockedRobustnessArg = microbialCompartmentBlockedRobustnessEnabled ? 'TRUE' : 'FALSE'
    """
set -euo pipefail
echo "infer_microbial_compartments.R md5: ${microbialCompartmentInferenceScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"

# Stability replicates are parallelized across R workers. Prevent each worker
# from creating an additional native numerical-library thread pool.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

Rscript "${microbialCompartmentInferenceScriptPath}" \\
  --prepared-dir "${prepared}" \\
  --outdir microbial_compartment_inference \\
  --sample-col "${microbialCompartmentSampleCol}" \\
  --k-min ${microbialCompartmentKMin} \\
  --k-max ${microbialCompartmentKMax} \\
  --min-cluster-size ${microbialCompartmentMinClusterSize} \\
  --min-cluster-fraction ${microbialCompartmentMinClusterFraction} \\
  --min-mean-silhouette ${microbialCompartmentMinMeanSilhouette} \\
  --min-stability-ari ${microbialCompartmentMinStabilityAri} \\
  --min-cluster-jaccard ${microbialCompartmentMinClusterJaccard} \\
  --stability-replicates ${microbialCompartmentStabilityReplicates} \\
  --stability-sample-fraction ${microbialCompartmentStabilitySampleFraction} \\
  --stability-primary-quantile ${microbialCompartmentStabilityPrimaryQuantile} \\
  --stability-near-tie-tolerance ${microbialCompartmentStabilityNearTieTolerance} \\
  --blocked-robustness-enabled ${blockedRobustnessArg} \\
  --prediction-strength-enabled ${predictionStrengthArg} \\
  --hierarchical-enabled ${hierarchicalArg} \\
  --seed ${microbialCompartmentSeed} \\
  --ncpus ${microbialCompartmentNcpus}

touch microbial_compartment_inference.done
"""
}

process MICROBIAL_COMPARTMENT_POSTHOC {
    cpus 1
    conda "${microbialCompartmentCondaEnvPath}"

    input:
    path(prepared)
    path(inference)
    path(metadata_table)
    path(taxonomy_table)

    output:
    path("microbial_compartment_posthoc_metadata.tsv"), emit: metadata
    path("microbial_compartment_posthoc.done"), emit: done

    when:
    microbialCompartmentEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
echo "posthoc_microbial_compartments.py md5: ${microbialCompartmentPosthocScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"
rm -rf microbial_compartment_staged

"\${CONDA_PREFIX}/bin/python" "${microbialCompartmentPosthocScriptPath}" \\
  --prepared-dir "${prepared}" \\
  --inference-dir "${inference}" \\
  --metadata "${metadata_table}" \\
  --taxonomy "${taxonomy_table}" \\
  --outdir microbial_compartment_staged \\
  --sample-col "${microbialCompartmentSampleCol}" \\
  --environmental-compartment-cols "${microbialCompartmentEnvironmentalCols}" \\
  --depth-col "${microbialCompartmentDepthCol}" \\
  --date-col "${microbialCompartmentDateCol}" \\
  --dominant-top-n ${microbialCompartmentDominantTopN}

"\${CONDA_PREFIX}/bin/python" -c 'from pathlib import Path; import shutil; src=Path("microbial_compartment_staged"); dst=Path("${microbialCompartmentOutputDirAbs}"); shutil.rmtree(dst, ignore_errors=True); dst.parent.mkdir(parents=True, exist_ok=True); shutil.copytree(src, dst)'
cp microbial_compartment_staged/tables/microbial_compartment_posthoc_metadata.tsv microbial_compartment_posthoc_metadata.tsv
touch microbial_compartment_posthoc.done
"""
}

process MICROBIAL_COMPARTMENT_PLOTS {
    cpus 1
    conda "${microbialCompartmentCondaEnvPath}"

    input:
    path(posthoc_done)

    output:
    path("microbial_compartment_plots.done"), emit: done

    when:
    microbialCompartmentEnabled && microbialCompartmentPlotEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
echo "plot_microbial_compartments.py md5: ${microbialCompartmentPlotScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"
mkdir -p "${microbialCompartmentOutputDirAbs}/plots"

"\${CONDA_PREFIX}/bin/python" "${microbialCompartmentPlotScriptPath}" \\
  --module-dir "${microbialCompartmentOutputDirAbs}" \\
  --outdir "${microbialCompartmentOutputDirAbs}/plots" \\
  --sample-col "${microbialCompartmentSampleCol}" \\
  --depth-col "${microbialCompartmentDepthCol}" \\
  --date-col "${microbialCompartmentDateCol}" \\
  --formats "${microbialCompartmentFormats}"

touch microbial_compartment_plots.done
"""
}

process MICROBIAL_STATE_INTERPRETATION {
    cpus 1
    conda "${microbialStateInterpretationCondaEnvPath}"

    input:
    path(microbial_compartment_done)
    path(modules_all)
    path(metadata_table)
    path(taxonomy_table)
    path(measurement_association_done)
    val(contour_time_step_days)
    val(contour_max_time_support_days)
    val(contour_max_depth_support_m)
    val(contour_reference_grid)

    output:
    path("microbial_state_interpretation.done"), emit: done

    when:
    microbialStateInterpretationEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def renewalArgs = microbialStateInterpretationRenewalEvents ? """  --renewal-events "${microbialStateInterpretationRenewalEvents}" \\
  --renewal-date-col "${microbialStateInterpretationRenewalDateCol}" \\
""" : ''
    def contourReferenceArg = contour_reference_grid ? """  --contour-reference-grid "${contour_reference_grid}" \\
""" : ''
    """
set -euo pipefail
echo "microbial_state_interpretation.py md5: ${microbialStateInterpretationScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"
rm -rf microbial_state_interpretation_staged

"\${CONDA_PREFIX}/bin/python" "${microbialStateInterpretationScriptPath}" \\
  --microbial-dir "${microbialCompartmentOutputDirAbs}" \\
  --metadata "${metadata_table}" \\
  --modules "${modules_all}" \\
  --taxonomy "${taxonomy_table}" \\
  --measurement-matrix "${measurementAssociationOutputDirAbs}/tables/measurement_matrix.samples_by_measurement.tsv" \\
  --outdir microbial_state_interpretation_staged \\
  --sample-col "${microbialStateInterpretationSampleCol}" \\
  --cruise-col "${microbialStateInterpretationCruiseCol}" \\
  --depth-col "${microbialStateInterpretationDepthCol}" \\
  --season-col "${microbialStateInterpretationSeasonCol}" \\
  --year-col "${microbialStateInterpretationYearCol}" \\
  --date-col "${microbialStateInterpretationDateCol}" \\
  --renewal-col "${microbialStateInterpretationRenewalCol}" \\
  --renewal-levels "${microbialStateInterpretationRenewalLevels}" \\
  --environmental-compartment-cols "${microbialStateInterpretationEnvironmentalCols}" \\
  --hybrid-col "${microbialStateInterpretationHybridCol}" \\
  --permutations ${microbialStateInterpretationPermutations} \\
  --bootstrap-replicates ${microbialStateInterpretationBootstrapReplicates} \\
  --seed ${microbialStateInterpretationSeed} \\
  --hybrid-palette "${microbialStateInterpretationHybridPalette}" \\
  --hybrid-order "${microbialStateInterpretationHybridOrder}" \\
  --mc-palette "${microbialStateInterpretationMcPalette}" \\
  --depth-palette "${microbialStateInterpretationDepthPalette}" \\
  --mc-order "${microbialStateInterpretationMcOrder}" \\
  --asv-top-n ${microbialStateInterpretationAsvTopN} \\
  --selected-asvs "${microbialStateInterpretationSelectedAsvs}" \\
  --asv-selection-metric "${microbialStateInterpretationAsvSelectionMetric}" \\
  --maximum-depth ${microbialStateInterpretationMaximumDepth} \\
  --depth-step ${microbialStateInterpretationDepthStep} \\
  --time-subdivisions-per-month ${microbialStateInterpretationTimeSubdivisions} \\
  --time-sigma-months ${microbialStateInterpretationTimeSigma} \\
  --depth-visual-sigma-m ${microbialStateInterpretationDepthSigma} \\
  --contour-time-step-days ${contour_time_step_days} \\
  --contour-max-time-support-days ${contour_max_time_support_days} \\
  --contour-max-depth-support-m ${contour_max_depth_support_m} \\
${contourReferenceArg}  --exclude-curtain-label-pattern "${microbialStateInterpretationExcludeCurtainPattern}" \\
${renewalArgs}  --formats "${microbialStateInterpretationFormats}"

"\${CONDA_PREFIX}/bin/python" -c 'from pathlib import Path; import shutil; src=Path("microbial_state_interpretation_staged"); dst=Path("${microbialStateInterpretationOutputDirAbs}"); shutil.rmtree(dst, ignore_errors=True); dst.parent.mkdir(parents=True, exist_ok=True); shutil.copytree(src, dst)'
[[ -f "${microbialStateInterpretationOutputDirAbs}/tables/microbial_state_curtain_renewal_onsets.tsv" ]] || { echo "Missing microbial-state curtain renewal-onset audit" >&2; exit 1; }
[[ -f "${microbialStateInterpretationOutputDirAbs}/tables/continuous_compartment_grid_microbial.tsv.gz" ]] || { echo "Missing BASIN-comparable microbial-compartment contour grid" >&2; exit 1; }
[[ -f "${microbialStateInterpretationOutputDirAbs}/plots/continuous_time_depth_compartments_microbial.pdf" ]] || { echo "Missing BASIN-comparable microbial-compartment contour PDF" >&2; exit 1; }
[[ -f "${microbialStateInterpretationOutputDirAbs}/tables/microbial_compartment_evidence_summary.tsv" ]] || { echo "Missing consolidated microbial-compartment evidence summary" >&2; exit 1; }
[[ -f "${microbialStateInterpretationOutputDirAbs}/tables/ecological_module_home_microbial_compartment.tsv" ]] || { echo "Missing ecological-module MC home table" >&2; exit 1; }
[[ -f "${microbialStateInterpretationOutputDirAbs}/tables/microbial_state_module_association_statistics.tsv" ]] || { echo "Missing depth-adjusted module association statistics" >&2; exit 1; }
[[ -f "${microbialStateInterpretationOutputDirAbs}/tables/microbial_compartment_renewal_association.tsv" ]] || { echo "Missing MC-renewal association statistics" >&2; exit 1; }
[[ -f "${microbialStateInterpretationOutputDirAbs}/tables/ecological_module_renewal_association.tsv" ]] || { echo "Missing module-renewal association statistics" >&2; exit 1; }
touch microbial_state_interpretation.done
"""
}

process COMMUNITY_TURNOVER_PREPARE {
    cpus 1
    conda "${communityTurnoverCondaEnvPath}"

    input:
    path(asv_counts)
    path(spieceasi_filter_audit)
    path(metadata_table)

    output:
    path("community_turnover_prepared"), emit: prepared
    path("community_turnover_prepare.done"), emit: done

    when:
    communityTurnoverEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def transposeArg = communityTurnoverTranspose ? '--transpose' : '--no-transpose'
    """
set -euo pipefail
echo "prepare_microbial_compartments.py md5: ${microbialCompartmentPrepareScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"

"\${CONDA_PREFIX}/bin/python" "${microbialCompartmentPrepareScriptPath}" \\
  --counts "${asv_counts}" \\
  --spieceasi-filter-audit "${spieceasi_filter_audit}" \\
  --metadata "${metadata_table}" \\
  --outdir community_turnover_prepared \\
  --sample-col "${communityTurnoverSampleCol}" \\
  --asv-id-col "${measurementAssociationAsvIdCol}" \\
  ${transposeArg} \\
  --minimum-prevalence ${communityTurnoverMinimumPrevalence} \\
  --minimum-max-relative-abundance ${communityTurnoverMinimumMaxRelativeAbundance} \\
  --zero-replacement "${communityTurnoverZeroReplacement}" \\
  --multiplicative-delta ${communityTurnoverMultiplicativeDelta} \\
  --pseudocount ${communityTurnoverPseudocount}

touch community_turnover_prepare.done
"""
}

process COMMUNITY_TURNOVER_ANALYSIS {
    cpus 1
    conda "${communityTurnoverCondaEnvPath}"

    input:
    path(prepared)
    path(metadata_table)

    output:
    path("community_turnover_staged"), emit: analysis
    path("community_turnover_analysis.done"), emit: done

    when:
    communityTurnoverEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
echo "community_turnover.py md5: ${communityTurnoverScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"

"\${CONDA_PREFIX}/bin/python" "${communityTurnoverScriptPath}" \\
  --prepared-dir "${prepared}" \\
  --metadata "${metadata_table}" \\
  --outdir community_turnover_staged \\
  --sample-col "${communityTurnoverSampleCol}" \\
  --profile-col "${communityTurnoverProfileCol}" \\
  --date-col "${communityTurnoverDateCol}" \\
  --depth-col "${communityTurnoverDepthCol}" \\
  --environmental-compartment-cols "${communityTurnoverEnvironmentalCols}" \\
  --primary-environmental-compartment-col "${communityTurnoverPrimaryEnvironmentalCol}" \\
  --distance-metrics "${communityTurnoverDistanceMetrics}" \\
  --fixed-depth-min-profile-fraction ${communityTurnoverFixedDepthMinProfileFraction} \\
  --minimum-time-difference-days ${communityTurnoverMinimumTimeDifferenceDays} \\
  --boundary-permutations ${communityTurnoverBoundaryPermutations} \\
  --boundary-bootstrap-replicates ${communityTurnoverBoundaryBootstrapReplicates} \\
  --seed ${communityTurnoverSeed}

touch community_turnover_analysis.done
"""
}

process COMMUNITY_TURNOVER_LCBD {
    cpus 1
    conda "${communityTurnoverCondaEnvPath}"

    input:
    path(prepared)
    path(metadata_table)
    path(turnover_analysis)

    output:
    path("community_turnover_lcbd.done"), emit: done

    when:
    communityTurnoverEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def withinProfileArg = communityTurnoverWithinProfileLcbdEnabled ? 'TRUE' : 'FALSE'
    """
set -euo pipefail
echo "lcbd_analysis.R md5: ${communityTurnoverLcbdScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"
rm -rf community_turnover_publish
cp -a "${turnover_analysis}" community_turnover_publish
cp "${prepared}/microbial_clustering_sample_audit.tsv" community_turnover_publish/audit/compositional_preparation_sample_audit.tsv
cp "${prepared}/microbial_clustering_asv_filter_audit.tsv" community_turnover_publish/audit/compositional_preparation_asv_filter_audit.tsv
cp "${prepared}/zero_replacement_audit.tsv" community_turnover_publish/audit/zero_replacement_audit.tsv
cp "${prepared}/microbial_clustering_preparation_config.json" community_turnover_publish/audit/compositional_preparation_config.json

Rscript "${communityTurnoverLcbdScriptPath}" \\
  --prepared-dir "${prepared}" \\
  --metadata "${metadata_table}" \\
  --outdir community_turnover_lcbd \\
  --sample-col "${communityTurnoverSampleCol}" \\
  --profile-col "${communityTurnoverProfileCol}" \\
  --date-col "${communityTurnoverDateCol}" \\
  --depth-col "${communityTurnoverDepthCol}" \\
  --environmental-compartment-cols "${communityTurnoverEnvironmentalCols}" \\
  --permutations ${communityTurnoverLcbdPermutations} \\
  --seed ${communityTurnoverSeed} \\
  --within-profile-enabled ${withinProfileArg} \\
  --within-profile-min-samples ${communityTurnoverWithinProfileLcbdMinSamples}

cp community_turnover_lcbd/lcbd.tsv community_turnover_publish/tables/lcbd.tsv
cp community_turnover_lcbd/lcbd_within_profile_audit.tsv community_turnover_publish/audit/lcbd_within_profile_audit.tsv
cp community_turnover_lcbd/lcbd_parameters.tsv community_turnover_publish/audit/lcbd_parameters.tsv
cp community_turnover_lcbd/lcbd_session_info.txt community_turnover_publish/audit/lcbd_session_info.txt
"\${CONDA_PREFIX}/bin/python" -c 'from pathlib import Path; import shutil; src=Path("community_turnover_publish"); dst=Path("${communityTurnoverOutputDirAbs}"); shutil.rmtree(dst, ignore_errors=True); dst.parent.mkdir(parents=True, exist_ok=True); shutil.copytree(src, dst)'
touch community_turnover_lcbd.done
"""
}

process COMMUNITY_TURNOVER_PLOTS {
    cpus 1
    conda "${communityTurnoverCondaEnvPath}"

    input:
    path(lcbd_done)

    output:
    path("community_turnover_plots.done"), emit: done

    when:
    communityTurnoverEnabled && communityTurnoverPlotEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
echo "plot_community_turnover.py md5: ${communityTurnoverPlotScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"
mkdir -p "${communityTurnoverOutputDirAbs}/plots"
"\${CONDA_PREFIX}/bin/python" "${communityTurnoverPlotScriptPath}" \\
  --module-dir "${communityTurnoverOutputDirAbs}" \\
  --outdir "${communityTurnoverOutputDirAbs}/plots" \\
  --primary-metric "${communityTurnoverPrimaryMetric}" \\
  --depth-col "${communityTurnoverDepthCol}" \\
  --date-col "${communityTurnoverDateCol}" \\
  --formats "${communityTurnoverFormats}"
touch community_turnover_plots.done
"""
}

process COMMUNITY_PREDICTOR_COMPARISON {
    cpus pipelineThreads
    conda "${communityPredictorCondaEnvPath}"

    input:
    path(metadata_table)
    path(asv_counts)
    val(community_predictor_script_hash)

    output:
    path("community_predictor_comparison.done"), emit: done
    path("sample_analysis_cohort.tsv"), optional: true, emit: sample_cohort

    when:
    communityPredictorEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
mkdir -p "${communityPredictorOutputDirAbs}"
export MPLCONFIGDIR="\$PWD/.mplconfig"
mkdir -p "\${MPLCONFIGDIR}"
echo "community_predictor_comparison.py md5: ${community_predictor_script_hash}"
# A rerun must never emit a stale cohort left in publication staging.
rm -f "${communityPredictorOutputDirAbs}/tables/sample_analysis_cohort.tsv"

"\${CONDA_PREFIX}/bin/python" "${communityPredictorScriptPath}" \
  --metadata "${metadata_table}" \
  --asv-counts "${asv_counts}" \
  --outdir "${communityPredictorOutputDirAbs}" \
  --sample-col "${communityPredictorSampleCol}" \
  --cruise-col "${communityPredictorCruiseCol}" \
  --year-col "${communityPredictorYearCol}" \
  --date-col "${communityPredictorDateCol}" \
  --season-col "${communityPredictorSeasonCol}" \
  --depth-col "${communityPredictorDepthCol}" \
  --cruise-depth-min-prevalence ${communityPredictorCruiseDepthMinPrevalence} \
  --pea-col "${communityPredictorPeaCol}" \
  --centroid-col "${communityPredictorCentroidCol}" \
  --cruise-group-col "${communityPredictorCruiseGroupCol}" \
  --cruise-group-probability-col "${communityPredictorCruiseGroupProbabilityCol}" \
  --cruise-group-uncertain-col "${communityPredictorCruiseGroupUncertainCol}" \
  --renewal-group-col "${communityPredictorRenewalGroupCol}" \
  --o2-group-col "${communityPredictorO2Col}" \
  --gmm-group-col "${communityPredictorGmmCol}" \
  --hybrid-group-col "${communityPredictorHybridCol}" \
  --assignment-source-col "${communityPredictorSourceCol}" \
  --observed-source-label "${communityPredictorObservedLabel}" \
  --min-group-n ${communityPredictorMinGroupN} \
  --permutations ${communityPredictorPermutations} \
  --seed ${communityPredictorSeed} \
  --formats "${communityPredictorFormats}"

# Emit a task-local copy so consumers depend on this run, not a published path.
if [[ -f "${communityPredictorOutputDirAbs}/tables/sample_analysis_cohort.tsv" ]]; then
  cp "${communityPredictorOutputDirAbs}/tables/sample_analysis_cohort.tsv" sample_analysis_cohort.tsv
fi
touch community_predictor_comparison.done
"""
}
process GROUPING_DIAGNOSTICS {
    cpus pipelineThreads
    conda "${groupingDiagnosticsCondaEnvPath}"

    input:
    path(metadata_table)
    path(asv_counts)

    output:
    path("grouping_diagnostics.done"), emit: done
    path("grouping_soft_label_assignments.tsv"), optional: true, emit: soft_assignments
    path("grouping_soft_label_validation.tsv"), optional: true, emit: soft_validation
    path("grouping_soft_label_validation_summary.tsv"), optional: true, emit: soft_validation_summary

    when:
    groupingDiagnosticsEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def groupColsArg = groupingDiagnosticsGroupCols.join(',')
    def cruiseGroupColsArg = groupingDiagnosticsCruiseGroupCols.join(',')
    def softLabelTargetColsArg = groupingDiagnosticsSoftLabelTargetCols.join(',')
    def softLabelExcludeArg = groupingDiagnosticsSoftLabelExcludeLabels.join(',')
    def softLabelArg = groupingDiagnosticsSoftLabelEnabled ? """  --soft-label-missing \\\n  --soft-label-k ${groupingDiagnosticsSoftLabelK} \\\n  --soft-label-group-cols "${softLabelTargetColsArg}" \\\n  --soft-label-exclude-labels "${softLabelExcludeArg}" \\\n  --soft-label-min-class-samples ${groupingDiagnosticsSoftLabelMinClassSamples} \\\n  --soft-label-distance-quantile ${groupingDiagnosticsSoftLabelDistanceQuantile} \\\n""" : ''
    def powerArg = groupingDiagnosticsPowerEnabled ? """  --power-enabled \\\n  --power-sample-sizes "${groupingDiagnosticsPowerSizes}" \\\n  --power-simulations ${groupingDiagnosticsPowerSimulations} \\\n  --power-permutations ${groupingDiagnosticsPowerPermutations} \\\n  --power-alpha ${groupingDiagnosticsPowerAlpha} \\\n""" : ''
"""
set -euo pipefail
mkdir -p "${groupingDiagnosticsOutputDirAbs}"
export MPLCONFIGDIR="\$PWD/.mplconfig"
mkdir -p "\${MPLCONFIGDIR}"
echo "grouping_diagnostics.py md5: ${groupingDiagnosticsScriptHash}"

"\${CONDA_PREFIX}/bin/python" "${groupingDiagnosticsScriptPath}" \\
  --metadata "${metadata_table}" \\
  --asv-counts "${asv_counts}" \\
  --outdir "${groupingDiagnosticsOutputDirAbs}" \\
  --sample-col "${groupingDiagnosticsSampleCol}" \\
  --group-cols "${groupColsArg}" \\
  --cruise-level-group-cols "${cruiseGroupColsArg}" \\
  --cruise-col "${groupingDiagnosticsCruiseCol}" \\
  --depth-col "${groupingDiagnosticsDepthCol}" \\
  --cruise-depth-min-prevalence ${groupingDiagnosticsCruiseDepthMinPrevalence} \\
  --baseline-group "${groupingDiagnosticsBaselineGroup}" \\
  --primary-group "${groupingDiagnosticsPrimaryGroup}" \\
  --group-palettes-json '${groupingDiagnosticsPaletteJson}' \\
  --group-orders-json '${groupingDiagnosticsOrderJson}' \\
  --metrics "${groupingDiagnosticsMetrics}" \\
  --transform "${groupingDiagnosticsTransform}" \\
  --permutations ${groupingDiagnosticsPermutations} \\
  --random-state ${groupingDiagnosticsRandomState} \\
  --formats "${groupingDiagnosticsFormats}" \\
${softLabelArg}${powerArg}  --power-min-groups ${groupingDiagnosticsPowerMinGroups}

if [[ -f "${groupingDiagnosticsOutputDirAbs}/tables/grouping_soft_label_assignments.tsv" ]]; then
  cp "${groupingDiagnosticsOutputDirAbs}/tables/grouping_soft_label_assignments.tsv" grouping_soft_label_assignments.tsv
  cp "${groupingDiagnosticsOutputDirAbs}/tables/grouping_soft_label_validation.tsv" grouping_soft_label_validation.tsv
  cp "${groupingDiagnosticsOutputDirAbs}/tables/grouping_soft_label_validation_summary.tsv" grouping_soft_label_validation_summary.tsv
fi

touch grouping_diagnostics.done
"""
}
process GROUP_LABEL_AUGMENTATION {
    cpus 1
    conda "${groupLabelAugmentationCondaEnvPath}"
    publishDir "${outputDir}/metadata", mode: 'copy', pattern: '*.augmented.tsv'

    input:
    path(metadata_table)
    path(asv_meta_table)
    path(soft_assignments)
    path(soft_validation_summary)

    output:
    path("metadata_updated_micro.augmented.tsv"), emit: metadata_augmented
    path("ASV_meta_micro.augmented.tsv"), emit: asv_meta_augmented
    path("group_label_augmentation_audit.tsv"), emit: audit
    path("group_label_augmentation.done"), emit: done

    when:
    groupingDiagnosticsApplySoftLabels

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def excludedLabelsArg = groupingDiagnosticsSoftLabelExcludeLabels.join(',')
    """
set -euo pipefail
echo "group_label_augmentation.py md5: ${groupLabelAugmentationScriptHash}"

"\${CONDA_PREFIX}/bin/python" "${groupLabelAugmentationScriptPath}" \\
  --metadata "${metadata_table}" \\
  --asv-meta "${asv_meta_table}" \\
  --assignments "${soft_assignments}" \\
  --validation-summary "${soft_validation_summary}" \\
  --sample-col "${groupingDiagnosticsSampleCol}" \\
  --target-col "${groupingDiagnosticsSoftLabelTargetCol}" \\
  --exclude-labels "${excludedLabelsArg}" \\
  --min-confidence ${groupingDiagnosticsSoftLabelMinConfidence} \\
  --min-neighbor-agreement ${groupingDiagnosticsSoftLabelMinNeighborAgreement} \\
  --min-cv-balanced-accuracy ${groupingDiagnosticsSoftLabelMinCvBalancedAccuracy}

mkdir -p "${groupingDiagnosticsOutputDirAbs}/tables"
cp group_label_augmentation_audit.tsv "${groupingDiagnosticsOutputDirAbs}/tables/group_label_augmentation_audit.tsv"
touch group_label_augmentation.done
"""
}
process GROUP_POWER_ANALYSIS {
    cpus pipelineThreads
    conda "${powerAnalysisCondaEnvPath}"

    input:
    path(asv_meta)
    path(asv_counts)
    val(indicspecies_ready)

    output:
    path("power_analysis.done"), emit: done

    when:
    powerAnalysisEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def skipEstimateFlag = powerAnalysisSkipEstimate ? '1' : '0'
    def skipPlotFlag = powerAnalysisSkipPlot ? '1' : '0'
    def keepContralateralFlag = powerAnalysisKeepContralateralInCancer ? '1' : '0'
    """
set -euo pipefail
mkdir -p "${powerAnalysisOutputDirAbs}"

SUMMARY_INPUT_DIR="power_analysis_inputs"
mkdir -p "\${SUMMARY_INPUT_DIR}"

"\${CONDA_PREFIX}/bin/python" "${masterSummaryScriptPath}" \\
  --data-dir "${outputDir}" \\
  --asv-meta "${asv_meta}" \\
  --asv-counts "${asv_counts}" \\
  --clustermaps-dir "${clustermapsOutputDirAbs}" \\
  --indicspecies-dir "${powerAnalysisIndicspeciesDir}" \\
  --spieceasi-dir "${spieceasiOutputDirAbs}" \\
  --outdir "\${SUMMARY_INPUT_DIR}" \\
  --max-direct-cols ${masterSummaryMaxDirectCols}

if [[ ! -f "\${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" ]]; then
  echo "Missing expected power-analysis input: \${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" >&2
  exit 1
fi

POWER_INDICSPECIES_ARGS=()
if [[ -d "${powerAnalysisIndicspeciesDir}" ]]; then
  POWER_INDICSPECIES_ARGS=(--indicspecies-dir "${powerAnalysisIndicspeciesDir}")
fi

POWER_SKIP_ARGS=()
if [[ "${skipEstimateFlag}" == "1" ]]; then
  POWER_SKIP_ARGS+=(--skip-estimate)
fi
if [[ "${skipPlotFlag}" == "1" ]]; then
  POWER_SKIP_ARGS+=(--skip-plot)
fi

POWER_CONTRALATERAL_ARGS=()
if [[ "${keepContralateralFlag}" == "1" ]]; then
  POWER_CONTRALATERAL_ARGS+=(--keep-contralateral-in-cancer)
fi

bash "${powerAnalysisScriptPath}" \\
  --data-long "\${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" \\
  --data-wide "${asv_counts}" \\
  --outdir "${powerAnalysisOutputDirAbs}" \\
  --sample-col "${powerAnalysisSampleCol}" \\
  --patient-col "${powerAnalysisPatientCol}" \\
  --case-col "${powerAnalysisCaseCol}" \\
  --type-col "${powerAnalysisTypeCol}" \\
  --sample-sizes-cancer "${powerAnalysisSampleSizesCancer}" \\
  --sample-sizes-stype "${powerAnalysisSampleSizesStype}" \\
  --n-simulations ${powerAnalysisNSimulations} \\
  --n-perm ${powerAnalysisNPerm} \\
  --alpha ${powerAnalysisAlpha} \\
  --seed ${powerAnalysisSeed} \\
  --transform "${powerAnalysisTransform}" \\
  --contralateral-sample-types "${powerAnalysisContralateralTypes}" \\
  "\${POWER_SKIP_ARGS[@]}" \\
  "\${POWER_CONTRALATERAL_ARGS[@]}" \\
  "\${POWER_INDICSPECIES_ARGS[@]}"

touch power_analysis.done
"""
}
process TAXONOMY_GROUP_ASSOCIATION {
    cpus pipelineThreads
    conda "${taxonomyPatientAwareCondaEnvPath}"

    input:
    path(asv_meta)
    path(asv_counts)

    output:
    path("taxonomy_patient_aware.done"), emit: done

    when:
    taxonomyPatientAwareEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def excludeContralateralFlag = taxonomyPatientAwareExcludeContralateral ? '1' : '0'
    def skipOmnibusFlag = taxonomyPatientAwareSkipOmnibus ? '1' : '0'
    def runComparisonFlag = taxonomyPatientAwareRunComparison ? '1' : '0'
    def cancerResultsArg = taxonomyPatientAwareRunComparison ?
        "--cancer-results \"${taxonomyPatientAwareOutputDirAbs}/cancer_vs_control/taxonomic_abundance_observed.tsv\"" :
        ''
    """
set -euo pipefail
mkdir -p "${taxonomyPatientAwareOutputDirAbs}"

SUMMARY_INPUT_DIR="taxonomy_patient_aware_inputs"
mkdir -p "\${SUMMARY_INPUT_DIR}"

"\${CONDA_PREFIX}/bin/python" "${masterSummaryScriptPath}" \\
  --data-dir "${outputDir}" \\
  --asv-meta "${asv_meta}" \\
  --asv-counts "${asv_counts}" \\
  --clustermaps-dir "${clustermapsOutputDirAbs}" \\
  --indicspecies-dir "${indicspeciesOutputDirAbs}" \\
  --spieceasi-dir "${spieceasiOutputDirAbs}" \\
  --outdir "\${SUMMARY_INPUT_DIR}" \\
  --max-direct-cols ${masterSummaryMaxDirectCols}

if [[ ! -f "\${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" ]]; then
  echo "Missing expected taxonomy patient-aware input: \${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" >&2
  exit 1
fi

TAXONOMY_CONTRALATERAL_ARGS=()
if [[ "${excludeContralateralFlag}" == "1" ]]; then
  TAXONOMY_CONTRALATERAL_ARGS+=(--exclude-contralateral-in-cancer)
else
  TAXONOMY_CONTRALATERAL_ARGS+=(--keep-contralateral-in-cancer)
fi

TAXONOMY_OMNIBUS_ARGS=()
if [[ "${skipOmnibusFlag}" == "1" ]]; then
  TAXONOMY_OMNIBUS_ARGS+=(--skip-omnibus)
fi

SAMPLETYPE_OUTDIR="${taxonomyPatientAwareOutputDirAbs}/sample_type"
FIGURES_OUTDIR="${taxonomyPatientAwareOutputDirAbs}/figures"
rm -rf "${taxonomyPatientAwareOutputDirAbs}/cancer_vs_control"
mkdir -p "\${SAMPLETYPE_OUTDIR}" "\${FIGURES_OUTDIR}"

if [[ "${runComparisonFlag}" == "1" ]]; then
  CANCER_OUTDIR="${taxonomyPatientAwareOutputDirAbs}/cancer_vs_control"
  mkdir -p "\${CANCER_OUTDIR}"
  "\${CONDA_PREFIX}/bin/python" "${taxonomicAbundanceObservedScriptPath}" \\
    --data-long "\${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" \\
    --tax-levels "${taxonomyPatientAwareTaxLevels}" \\
    --sample-types "${taxonomyPatientAwareSampleTypes}" \\
    --case-groups "${taxonomyPatientAwareComparisonGroups}" \\
    --sample-col "${taxonomyPatientAwareSampleCol}" \\
    --patient-col "${taxonomyPatientAwarePatientCol}" \\
    --case-col "${taxonomyPatientAwareCaseCol}" \\
    --type-col "${taxonomyPatientAwareTypeCol}" \\
    --count-col "${taxonomyPatientAwareCountCol}" \\
    --min-prevalence ${taxonomyPatientAwareMinPrevalence} \\
    --contralateral-col "${taxonomyPatientAwareContralateralCol}" \\
    --cancer-site-col "${taxonomyPatientAwareCancerSiteCol}" \\
    --lung-side-col "${taxonomyPatientAwareLungSideCol}" \\
    --contralateral-value "${taxonomyPatientAwareContralateralValue}" \\
    --contralateral-sample-types "${taxonomyPatientAwareContralateralTypes}" \\
    --transform "${taxonomyPatientAwareTransform}" \\
    --outdir "\${CANCER_OUTDIR}" \\
    "\${TAXONOMY_CONTRALATERAL_ARGS[@]}"
fi

"\${CONDA_PREFIX}/bin/python" "${taxonomicSampleTypeObservedScriptPath}" \\
  --data-long "\${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" \\
  --tax-levels "${taxonomyPatientAwareTaxLevels}" \\
  --sample-types "${taxonomyPatientAwareSampleTypes}" \\
  --sample-col "${taxonomyPatientAwareSampleCol}" \\
  --patient-col "${taxonomyPatientAwarePatientCol}" \\
  --case-col "${taxonomyPatientAwareCaseCol}" \\
  --type-col "${taxonomyPatientAwareTypeCol}" \\
  --count-col "${taxonomyPatientAwareCountCol}" \\
  --min-prevalence ${taxonomyPatientAwareMinPrevalence} \\
  --contralateral-col "${taxonomyPatientAwareContralateralCol}" \\
  --cancer-site-col "${taxonomyPatientAwareCancerSiteCol}" \\
  --lung-side-col "${taxonomyPatientAwareLungSideCol}" \\
  --contralateral-value "${taxonomyPatientAwareContralateralValue}" \\
  --contralateral-sample-types "${taxonomyPatientAwareContralateralTypes}" \\
  --transform "${taxonomyPatientAwareTransform}" \\
  --outdir "\${SAMPLETYPE_OUTDIR}" \\
  "\${TAXONOMY_CONTRALATERAL_ARGS[@]}" \\
  "\${TAXONOMY_OMNIBUS_ARGS[@]}"

"\${CONDA_PREFIX}/bin/python" "${plotTaxonomicObservedScriptPath}" \\
  --data-long "\${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" \\
  ${cancerResultsArg} \\
  --sampletype-results "\${SAMPLETYPE_OUTDIR}/taxonomic_sample_type_observed_pairwise.tsv" \\
  --outdir "\${FIGURES_OUTDIR}" \\
  --alpha ${taxonomyPatientAwareAlpha} \\
  --top-n ${taxonomyPatientAwareTopN} \\
  --sample-col "${taxonomyPatientAwareSampleCol}" \\
  --patient-col "${taxonomyPatientAwarePatientCol}" \\
  --type-col "${taxonomyPatientAwareTypeCol}" \\
  --case-col "${taxonomyPatientAwareCaseCol}" \\
  --count-col "${taxonomyPatientAwareCountCol}" \\
  --type-palette "${taxonomyPatientAwareTypePalette}" \\
  --case-palette "${taxonomyPatientAwareCasePalette}"

touch taxonomy_patient_aware.done
"""
}
process PAIRED_GROUP_CONTRAST {
    cpus pipelineThreads
    conda "${lungStatusAnalysisCondaEnvPath}"

    input:
    path(asv_meta)
    path(asv_counts)

    output:
    path("lung_status_analysis.done"), emit: done

    when:
    lungStatusAnalysisEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
mkdir -p "${lungStatusAnalysisOutputDirAbs}"

SUMMARY_INPUT_DIR="lung_status_analysis_inputs"
mkdir -p "\${SUMMARY_INPUT_DIR}"

"\${CONDA_PREFIX}/bin/python" "${masterSummaryScriptPath}" \\
  --data-dir "${outputDir}" \\
  --asv-meta "${asv_meta}" \\
  --asv-counts "${asv_counts}" \\
  --clustermaps-dir "${clustermapsOutputDirAbs}" \\
  --indicspecies-dir "${indicspeciesOutputDirAbs}" \\
  --spieceasi-dir "${spieceasiOutputDirAbs}" \\
  --outdir "\${SUMMARY_INPUT_DIR}" \\
  --max-direct-cols ${masterSummaryMaxDirectCols}

if [[ ! -f "\${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" ]]; then
  echo "Missing expected lung-status input: \${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" >&2
  exit 1
fi

IFS=',' read -r -a LUNG_SAMPLE_TYPES <<< "${lungStatusAnalysisSampleTypes}"
for sample_type in "\${LUNG_SAMPLE_TYPES[@]}"; do
  sample_type="\$(printf '%s' "\${sample_type}" | xargs)"
  [[ -n "\${sample_type}" ]] || continue
  sample_slug="\${sample_type// /_}"
  sample_root="${lungStatusAnalysisOutputDirAbs}/\${sample_slug}"
  data_dir="\${sample_root}/data"
  results_dir="\${sample_root}/results"
  figures_dir="\${sample_root}/figures"
  mkdir -p "\${data_dir}" "\${results_dir}" "\${figures_dir}"

  "\${CONDA_PREFIX}/bin/python" "${prepareLungStatusScriptPath}" \\
    --input "\${SUMMARY_INPUT_DIR}/ASV_master_long.tsv" \\
    --sample-type "\${sample_type}" \\
    --sample-col "${lungStatusAnalysisSampleCol}" \\
    --type-col "${lungStatusAnalysisTypeCol}" \\
    --case-col "${lungStatusAnalysisCaseCol}" \\
    --patient-col "${lungStatusAnalysisPatientCol}" \\
    --cancer-site-col "${lungStatusAnalysisCancerSiteCol}" \\
    --lung-code-col "${lungStatusAnalysisLungCodeCol}" \\
    --tumor-side-col "${lungStatusAnalysisTumorSideCol}" \\
    --contralateral-col "${lungStatusAnalysisContralateralCol}" \\
    --healthy-col "${lungStatusAnalysisHealthyCol}" \\
    --lung-status-col "${lungStatusAnalysisStatusCol}" \\
    --status-a-value "${lungStatusAnalysisStatusAValue}" \\
    --status-b-value "${lungStatusAnalysisStatusBValue}" \\
    --reference-status-value "${lungStatusAnalysisReferenceStatusValue}" \\
    --outdir "\${data_dir}"

  meta_file="\$(find "\${data_dir}" -maxdepth 1 -name '*_metadata.tsv' | head -n 1)"
  asv_file="\$(find "\${data_dir}" -maxdepth 1 -name '*_ASV_table.tsv' | head -n 1)"
  if [[ -z "\${meta_file}" || -z "\${asv_file}" ]]; then
    echo "Missing prepared lung-status inputs for sample type: \${sample_type}" >&2
    exit 1
  fi

  Rscript "${lungStatusAnalysisScriptPath}" "\${meta_file}" "\${asv_file}" "\${results_dir}" \
    ${lungStatusAnalysisPermutations} ${lungStatusAnalysisSeed}

  "\${CONDA_PREFIX}/bin/python" "${plotLungStatusScriptPath}" \\
    --metadata "\${meta_file}" \\
    --patient-level "\${results_dir}/patient_level_metadata.tsv" \\
    --distances "\${results_dir}/patient_level_bray_distances.tsv" \\
    --summary "\${results_dir}/lung_status_contrasts_summary.tsv" \\
    --pairdist-a "\${results_dir}/contrast_A_pairwise_distances.tsv" \\
    --asv-table "\${asv_file}" \\
    --outdir "\${figures_dir}"
done

touch lung_status_analysis.done
"""
}
