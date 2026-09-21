process CLUSTERMAPS {
    cpus pipelineThreads
    conda "${clustermapsCondaEnvPath}"

    input:
    path(asv_meta)
    path(metadata_table)
    path(isa_table)

    output:
    path("clustermaps.done"), emit: done

    when:
    clustermapsEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def group3ColArg = clustermapsGroup3Col ? """  --group3-col "${clustermapsGroup3Col}" \\\n""" : ''
    def excludeGroup1Arg = clustermapsExcludeGroup1 ? """  --exclude-group1 "${clustermapsExcludeGroup1}" \\\n""" : ''
    def group1PaletteArg = clustermapsGroup1Palette ? """  --group1-palette "${clustermapsGroup1Palette}" \\\n""" : ''
    def group2PaletteArg = clustermapsGroup2Palette ? """  --group2-palette "${clustermapsGroup2Palette}" \\\n""" : ''
    def group3PaletteArg = clustermapsGroup3Palette ? """  --group3-palette "${clustermapsGroup3Palette}" \\\n""" : ''
    def clustermapsGroup1OrderArg = clustermapsGroup1Order && !clustermapsGroup1Order.isEmpty() ? """  --group1-order "${clustermapsGroup1Order.join(',')}" \\\n""" : ''
    def clustermapsRunMitoFlag = clustermapsRunMito ? '1' : '0'
    def isaTableArg = isa_table.name == 'empty_modules.tsv' ? '' : isa_table.toString()
    def isaSigColsArg = clustermapsIsaSignificanceCols ? """  --isa-significance-cols "${clustermapsIsaSignificanceCols}" \\\n""" : ''
    def isaStatColsArg = clustermapsIsaStatCols ? """  --isa-stat-cols "${clustermapsIsaStatCols}" \\\n""" : ''
    def clustermapsFormatsArg = clustermapsFormats ? """  --formats "${clustermapsFormats}" \\\n""" : ''
    def clustermapsFigWidthArg = clustermapsFigWidth != null ? """  --figwidth ${clustermapsFigWidth} \\\n""" : ''
    def clustermapsRowHeightArg = clustermapsRowHeight != null ? """  --row-height ${clustermapsRowHeight} \\\n""" : ''
    def clustermapsMinHeightArg = clustermapsMinHeight != null ? """  --min-height ${clustermapsMinHeight} \\\n""" : ''
    def clustermapsMaxHeightArg = clustermapsMaxHeight != null ? """  --max-height ${clustermapsMaxHeight} \\\n""" : ''
    def clustermapsSkipAsvPlotArg = clustermapsPlotAsvLevel ? '' : "  --skip-asv-plot \\\n"
    """
set -euo pipefail
mkdir -p "${clustermapsOutputDirAbs}"
mkdir -p "${clustermapsMitoOutputDirAbs}"

ISA_ARGS=()
if [[ -n "${isaTableArg}" ]]; then
  if [[ ! -f "${isaTableArg}" ]]; then
    echo "[STOP] Staged CLUSTERMAPS ISA input is missing: ${isaTableArg}" >&2
    exit 1
  fi
  echo "[i] CLUSTERMAPS using staged ISA table: ${isaTableArg}" >&2
  ISA_ARGS=(--isa "${isaTableArg}")
fi

if [[ "${clustermapsRunMitoFlag}" == "1" && -f "${clustermapsMitoInputPath}" ]]; then
  "\${CONDA_PREFIX}/bin/python" "${clustermapsScriptPath}" \\
    --asv-meta "${asv_meta}" \\
    --metadata "${metadata_table}" \\
    --outdir "${clustermapsOutputDirAbs}" \\
    --sample-col "${clustermapsSampleCol}" \\
    --sample-code-col "${clustermapsSampleCodeCol}" \\
    --asv-id-col "${clustermapsAsvIdCol}" \\
    --group1-col "${clustermapsGroup1Col}" \\
    --group2-col "${clustermapsGroup2Col}" \\
${group3ColArg}${clustermapsGroup1OrderArg}${excludeGroup1Arg}${group1PaletteArg}${group2PaletteArg}\
${group3PaletteArg}    --ranks "${clustermapsRanks}" \\
    --topN "${clustermapsTopN}" \\
    --count-col "${clustermapsCountCol}" \\
    --isa-min-stat ${clustermapsIsaMinStat} \\
${clustermapsFormatsArg}${clustermapsFigWidthArg}${clustermapsRowHeightArg}${clustermapsMinHeightArg}${clustermapsMaxHeightArg}\
${isaSigColsArg}${isaStatColsArg}${clustermapsSkipAsvPlotArg}    --mito-sample-mode "${clustermapsMitoSampleMode}" \\
    --mito-asv "${clustermapsMitoInputPath}" \\
    --mito-outdir "${clustermapsMitoOutputDirAbs}" \\
    \${ISA_ARGS[@]}
else
  "\${CONDA_PREFIX}/bin/python" "${clustermapsScriptPath}" \\
    --asv-meta "${asv_meta}" \\
    --metadata "${metadata_table}" \\
    --outdir "${clustermapsOutputDirAbs}" \\
    --sample-col "${clustermapsSampleCol}" \\
    --sample-code-col "${clustermapsSampleCodeCol}" \\
    --asv-id-col "${clustermapsAsvIdCol}" \\
    --group1-col "${clustermapsGroup1Col}" \\
    --group2-col "${clustermapsGroup2Col}" \\
${group3ColArg}${clustermapsGroup1OrderArg}${excludeGroup1Arg}${group1PaletteArg}${group2PaletteArg}\
${group3PaletteArg}    --ranks "${clustermapsRanks}" \\
    --topN "${clustermapsTopN}" \\
    --count-col "${clustermapsCountCol}" \\
    --isa-min-stat ${clustermapsIsaMinStat} \\
${clustermapsFormatsArg}${clustermapsFigWidthArg}${clustermapsRowHeightArg}${clustermapsMinHeightArg}${clustermapsMaxHeightArg}\
${isaSigColsArg}${isaStatColsArg}${clustermapsSkipAsvPlotArg}    --mito-sample-mode "${clustermapsMitoSampleMode}" \\
    \${ISA_ARGS[@]}
fi

touch clustermaps.done
"""
}
process SPIECEASI {
    cpus pipelineThreads
    conda "${spieceasiCondaEnvPath}"

    input:
    path(asv_counts)
    path(force_keep_asvs)
    path(asv_mag_links)

    output:
    path("spieceasi_network_pos_all.graphml"), emit: graph_all
    path("spieceasi_network_pos_thr.graphml"), emit: graph_thr
    path("spieceasi_node_features.csv"), emit: node_features
    path("spieceasi_filtering_audit.csv"), emit: filter_audit
    path("spieceasi.done"), emit: done

    when:
    spieceasiEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def transposeFlag = spieceasiTranspose ? 'TRUE' : 'FALSE'
    def removeZeroVarFlag = spieceasiRemoveZeroVar ? 'TRUE' : 'FALSE'
    def keepNegativeFlag = spieceasiKeepNegative ? 'TRUE' : 'FALSE'
    def forceFilterFlag = spieceasiForceFilter ? 'TRUE' : 'FALSE'
    def forceSpieceasiFlag = spieceasiForceSpieceasi ? 'TRUE' : 'FALSE'
    def forceGraphsFlag = spieceasiForceGraphs ? 'TRUE' : 'FALSE'
    def forceKeepAsvsArg = (indicspeciesEnabled && spieceasiForceKeepIsaAsvs) ? """  --force-keep-asvs "${force_keep_asvs}" \\\n""" : ''
    def forceKeepAsvMagArg = (asvMagLinkEnabled && spieceasiForceKeepAsvMagAsvs) ? """  --force-keep-asv-mag-links "${asv_mag_links}" \\
  --asv-mag-min-pident ${asvMagNetworkMinPident} \\
  --asv-mag-min-qcov ${asvMagNetworkMinQcov} \\
""" : ''
    """
set -euo pipefail
mkdir -p "${spieceasiOutputDirAbs}"

# pulsar parallelizes across R worker processes. Keep numerical libraries
# single-threaded inside each worker to avoid nested CPU oversubscription.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export RCPP_PARALLEL_NUM_THREADS=1

Rscript "${spieceasiScriptPath}" \\
  --counts "${asv_counts}" \\
  --outdir "${spieceasiOutputDirAbs}" \\
  --prefix "${spieceasiPrefix}" \\
  --transpose ${transposeFlag} \\
  --min-rel-abund ${spieceasiMinRelAbund} \\
  --min-prevalence ${spieceasiMinPrevalence} \\
${forceKeepAsvsArg}${forceKeepAsvMagArg}  --remove-zero-var ${removeZeroVarFlag} \\
  --method "${spieceasiMethod}" \\
  --lambda-min-ratio ${spieceasiLambdaMinRatio} \\
  --nlambda ${spieceasiNlambda} \\
  --rep-num ${spieceasiRepNum} \\
  --thresh ${spieceasiThresh} \\
  --pulsar-criterion "${spieceasiPulsarCriterion}" \\
  --ncores ${spieceasiNcores} \\
  --seed ${spieceasiSeed} \\
  --edge-threshold ${spieceasiEdgeThreshold} \\
  --keep-negative ${keepNegativeFlag} \\
  --layout-iters ${spieceasiLayoutIters} \\
  --force-filter ${forceFilterFlag} \\
  --force-spieceasi ${forceSpieceasiFlag} \\
  --force-graphs ${forceGraphsFlag}

GRAPH_ALL="${spieceasiOutputDirAbs}/${spieceasiPrefix}_network_pos_all.graphml"
GRAPH_THR="${spieceasiOutputDirAbs}/${spieceasiPrefix}_network_pos_thr.graphml"
NODE_FEATURES="${spieceasiOutputDirAbs}/${spieceasiPrefix}_node_features.csv"
FILTER_AUDIT="${spieceasiOutputDirAbs}/${spieceasiPrefix}_filtering_audit.csv"

for f in "\${GRAPH_ALL}" "\${GRAPH_THR}" "\${NODE_FEATURES}" "\${FILTER_AUDIT}"; do
  if [[ ! -f "\${f}" ]]; then
    echo "Missing expected SPIEC-EASI output: \${f}" >&2
    exit 1
  fi
done

ln -sf "\${GRAPH_ALL}" spieceasi_network_pos_all.graphml
ln -sf "\${GRAPH_THR}" spieceasi_network_pos_thr.graphml
ln -sf "\${NODE_FEATURES}" spieceasi_node_features.csv
ln -sf "\${FILTER_AUDIT}" spieceasi_filtering_audit.csv
touch spieceasi.done
"""
}
process NETWORK_MODULES {
    cpus pipelineThreads
    conda "${networkModulesCondaEnvPath}"

    input:
    path(graph_all, stageAs: 'network_graph_all.graphml')
    path(graph_thr, stageAs: 'network_graph_sub.graphml')

    output:
    path("network_modules_sub.tsv"), emit: modules_sub
    path("network_modules_all.tsv"), emit: modules_all
    path("network_modules_summary.tsv"), emit: summary
    path("network_modules_runs.tsv"), emit: runs
    path("network_modules.done"), emit: done

    when:
    networkEnabled && networkModulesEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def methodsCsv = networkModuleMethods.join(',')
    def resolutionsCsv = networkModuleResolutions.join(',')
    """
set -euo pipefail
mkdir -p "${spieceasiOutputDirAbs}"

Rscript "${networkModulesScriptPath}" \\
  --graph-sub "${graph_thr}" \\
  --graph-all "${graph_all}" \\
  --outdir "${spieceasiOutputDirAbs}" \\
  --prefix "${spieceasiPrefix}" \\
  --methods "${methodsCsv}" \\
  --primary-method "${networkModulePrimaryMethod}" \\
  --reps ${networkModuleReps} \\
  --resolutions "${resolutionsCsv}" \\
  --consensus-threshold ${networkModuleConsensusThreshold} \\
  --seed ${networkModuleSeed}

MODULES_SUB="${spieceasiOutputDirAbs}/${spieceasiPrefix}_modules_sub.tsv"
MODULES_ALL="${spieceasiOutputDirAbs}/${spieceasiPrefix}_modules_all.tsv"
MODULE_SUMMARY="${spieceasiOutputDirAbs}/${spieceasiPrefix}_module_summary.tsv"
MODULE_RUNS="${spieceasiOutputDirAbs}/${spieceasiPrefix}_module_runs.tsv"

for f in "\${MODULES_SUB}" "\${MODULES_ALL}" "\${MODULE_SUMMARY}" "\${MODULE_RUNS}"; do
  if [[ ! -f "\${f}" ]]; then
    echo "Missing expected network module output: \${f}" >&2
    exit 1
  fi
done

ln -sf "\${MODULES_SUB}" network_modules_sub.tsv
ln -sf "\${MODULES_ALL}" network_modules_all.tsv
ln -sf "\${MODULE_SUMMARY}" network_modules_summary.tsv
ln -sf "\${MODULE_RUNS}" network_modules_runs.tsv
touch network_modules.done
"""
}

process GENOME_COOCCURRENCE {
    cpus pipelineThreads
    conda "${genomeCooccurrenceCondaEnvPath}"

    input:
    path(genome_qc)
    path(metagenome_abundance, stageAs: 'metagenome_abundance_long.tsv')
    path(metatranscriptome_abundance, stageAs: 'metatranscriptome_abundance_long.tsv')
    path(sample_metadata, stageAs: 'genome_cooccurrence_metadata.tsv')
    path(asv_mag_mappings)
    path(asv_modules)
    path(isa_tables)
    path(titan_taxon_results)
    path(microbial_state_done)

    output:
    path("genome_cooccurrence.done"), emit: done

    when:
    genomeCooccurrenceEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
echo "genome proportionality scripts md5: ${genomeCooccurrenceScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"
rm -rf genome_cooccurrence_staged
mkdir -p genome_cooccurrence_staged/prepared genome_cooccurrence_staged/metagenome genome_cooccurrence_staged/metatranscriptome

"\${CONDA_PREFIX}/bin/python" "${genomeCooccurrencePrepareScriptPath}" \
  --genome-qc "${genome_qc}" \
  --metagenome-abundance "${metagenome_abundance}" \
  --metatranscriptome-abundance "${metatranscriptome_abundance}" \
  --asv-mag-mappings "${asv_mag_mappings}" \
  --asv-modules "${asv_modules}" \
  --outdir genome_cooccurrence_staged/prepared \
  --sample-col "${genomeCooccurrenceSampleCol}" \
  --genome-col "${genomeCooccurrenceGenomeCol}" \
  --metagenome-value-col "${genomeCooccurrenceMetagenomeValueCol}" \
  --metatranscriptome-value-col "${genomeCooccurrenceMetatranscriptomeValueCol}" \
  --input-scale "${genomeCooccurrenceInputScale}" \
  --input-feature-level "${genomeCooccurrenceInputFeatureLevel}" \
  --expected-species ${genomeCooccurrenceExpectedSpecies} \
  --minimum-read-count ${genomeCooccurrenceMinimumReadCount} \
  --closure-total ${genomeCooccurrenceClosureTotal}

for modality in metagenome metatranscriptome; do
  prefix="genome_\${modality}"
  outdir="genome_cooccurrence_staged/\${modality}"
  "\${CONDA_PREFIX}/bin/python" "${genomeCooccurrenceInferenceScriptPath}" \
    --matrix "genome_cooccurrence_staged/prepared/\${modality}_species_matrix.tsv" \
    --metadata "${sample_metadata}" \
    --outdir "\${outdir}" \
    --prefix "\${prefix}" \
    --metadata-sample-col "${genomeCooccurrenceMetadataSampleCol}" \
    --cruise-col "${genomeCooccurrenceCruiseCol}" \
    --season-col "${genomeCooccurrenceSeasonCol}" \
    --depth-col "${genomeCooccurrenceDepthCol}" \
    --month-col "${genomeCooccurrenceMonthCol}" \
    --sample-id-regex '${genomeCooccurrenceSampleIdRegex}' \
    --zero-replacement-fraction ${genomeCooccurrenceZeroReplacementFraction} \
    --min-relative-abundance ${genomeCooccurrenceMinRelAbund} \
    --min-prevalence ${genomeCooccurrenceMinPrevalence} \
    --min-abs-rho ${genomeCooccurrenceMinAbsRho} \
    --bootstrap-iterations ${genomeCooccurrenceBootstrapIterations} \
    --min-bootstrap-recovery ${genomeCooccurrenceMinBootstrapRecovery} \
    --min-sign-consistency ${genomeCooccurrenceMinSignConsistency} \
    --permutations ${genomeCooccurrencePermutations} \
    --permutation-strata "${genomeCooccurrencePermutationStrata}" \
    --max-q ${genomeCooccurrenceMaxQ} \
    --seed ${genomeCooccurrenceSeed}

  "\${CONDA_PREFIX}/bin/python" "${genomeCooccurrencePlotScriptPath}" \
    --graph "\${outdir}/\${prefix}_network_all.graphml" \
    --node-features "\${outdir}/\${prefix}_node_features.csv" \
    --network-clusters "\${outdir}/\${prefix}_proportionality_clusters.tsv" \
    --matrix "genome_cooccurrence_staged/prepared/\${modality}_species_matrix.tsv" \
    --pair-statistics "\${outdir}/\${prefix}_all_pair_statistics.tsv" \
    --node-metadata genome_cooccurrence_staged/prepared/representative_node_metadata.tsv \
    --asv-crosswalk genome_cooccurrence_staged/prepared/representative_asv_crosswalk.tsv \
    --hybrid-isa o2_subcompartment_final_indicator_species_summary.tsv \
    --mc-isa microbial_compartment_indicator_species_summary.tsv \
    --titan "${titan_taxon_results}" \
    --module-renewal-association "${microbialStateInterpretationOutputDirAbs}/tables/ecological_module_renewal_association.tsv" \
    --module-renewal-profiles "${microbialStateInterpretationOutputDirAbs}/tables/ecological_module_renewal_profiles.tsv" \
    --renewal-max-q 0.05 \
    --modality "\${modality}" \
    --outdir "\${outdir}/overlays" \
    --seed ${genomeCooccurrenceSeed}
done

"\${CONDA_PREFIX}/bin/python" -c 'from pathlib import Path; import shutil; src=Path("genome_cooccurrence_staged"); dst=Path("${genomeCooccurrenceOutputDirAbs}"); shutil.rmtree(dst, ignore_errors=True); dst.parent.mkdir(parents=True, exist_ok=True); shutil.copytree(src, dst)'

for modality in metagenome metatranscriptome; do
  prefix="genome_\${modality}"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/genome_network_node_annotations.tsv"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/genome_network_summary.tsv"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/network_phylum.pdf"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/network_propagated_asv_ecological_module.pdf"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/network_propagated_microbial_compartment.pdf"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/network_propagated_renewal_associated_modules.pdf"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/network_propagated_renewal_post_renewal_vs_stagnation.pdf"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/network_proportionality_clusters.pdf"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/proportionality_heatmap_all_pairs.pdf"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/overlays/proportionality_heatmap_supported_edges.pdf"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/\${prefix}_all_pair_statistics.tsv"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/\${prefix}_inference_summary.tsv"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/\${prefix}_software_versions.tsv"
  test -s "${genomeCooccurrenceOutputDirAbs}/\${modality}/\${prefix}_proportionality_clusters.tsv"
done
touch genome_cooccurrence.done
"""
}
process GRAPH_NETWORK {
    cpus pipelineThreads
    conda "${networkCondaEnvPath}"

    input:
    path(graph_all, stageAs: 'network_graph_all.graphml')
    path(graph_thr, stageAs: 'network_graph_sub.graphml')
    path(node_features)
    path(asv_counts)
    path(metadata_table)
    path(dep_asv_mag, stageAs: 'dep_asv_mag.tsv')
    path(taxonomy_table)
    path(indicspecies_tables)
    path(modules_sub, stageAs: 'network_modules_sub.tsv')
    path(modules_all, stageAs: 'network_modules_all.tsv')
    path(microbial_state_done, stageAs: 'microbial_state_interpretation.done')

    output:
    path("network.done"), emit: done
    path("network_modules_best_stats_all.tsv"), emit: module_best_stats_all
    path("spieceasi_network_layout_all.tsv"), emit: layout_all

    when:
    networkEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def networkModesArg = networkModes && !networkModes.isEmpty() ? """  --modes ${networkModes.collect { "\"${it}\"" }.join(' ')} \\\n""" : ''
    def networkModuleBestOnlyArg = networkModuleBestOnly ? """  --module-best-only \\\n""" : ''
    def networkModuleIsaOnlyArg = networkModuleIsaOnly ? """  --module-isa-only \\\n""" : ''
    def networkModuleColorByIsaArg = networkModuleColorByIsa ? """  --module-color-by-isa \\\n""" : ''
    def asvMagPairingArg = asvMagLinkEnabled ? """  --asv-mag-pairing "${asvMagLinkOutputDirAbs}/tables/asv2mag_pairing.tsv" \\\n""" : ''
    def moduleSubnetworksArg = networkModuleSubnetworksEnabled ? """  --module-subnetworks \\\n  --anchor-top-n ${asvMagNetworkAnchorTopN} \\\n  --module-subnetwork-prominence-metrics "${networkModuleSubnetworkProminenceMetrics}" \\\n  --module-subnetwork-prominence-threshold ${networkModuleSubnetworkProminenceThreshold} \\\n  --module-subnetwork-label-top-n ${networkModuleSubnetworkLabelTopN} \\\n  --module-subnetwork-prominence-min-area ${networkModuleSubnetworkProminenceMinArea} \\\n  --module-subnetwork-prominence-max-area ${networkModuleSubnetworkProminenceMaxArea} \\\n""" : ''
    def renewalAssociationArg = microbialStateInterpretationEnabled ? """  --module-renewal-association "${microbialStateInterpretationOutputDirAbs}/tables/ecological_module_renewal_association.tsv" \\
  --module-renewal-profiles "${microbialStateInterpretationOutputDirAbs}/tables/ecological_module_renewal_profiles.tsv" \\
  --module-renewal-max-q 0.05 \\
""" : ''
    def isaSummaryModeArg = 'default'
    """
set -euo pipefail
mkdir -p "${spieceasiOutputDirAbs}"

"\${CONDA_PREFIX}/bin/python" "${graphNetworkScriptPath}" \\
  --data-dir "${outputDir}" \\
  --outdir "${spieceasiOutputDirAbs}" \\
  --graph-pos-all "${graph_all}" \\
  --graph-pos-sub "${graph_thr}" \\
  --node-features "${node_features}" \\
  --asv-counts "${asv_counts}" \\
  --taxonomy "${taxonomy_table}" \\
  --metadata "${metadata_table}" \\
  --sample-col "${indicspeciesSampleCol}" \\
  --isa-group-cols "${networkIsaOverlayGroupsCsv}" \\
  --isa-summary-mode "${isaSummaryModeArg}" \\
  --isa-palette-map-json '${networkGroupPaletteJson}' \\
  --isa-order-map-json '${networkGroupOrderJson}' \\
  --isa-focus-map-json '${networkFocusLabelJson}' \\
${asvMagPairingArg}\
  --color-col "${networkColorCol}" \\
${networkModuleBestOnlyArg}${networkModuleIsaOnlyArg}${networkModuleColorByIsaArg}  --module-best-min-size ${networkModuleBestMinSize} \\
  --module-best-min-stability ${networkModuleBestMinStability} \\
  --module-best-top-n ${networkModuleBestTopN} \\
  --module-isa-source "${networkModuleIsaSource}" \\
  --module-isa-min-stat ${networkModuleIsaMinStat} \\
  --module-isa-max-q ${networkModuleIsaMaxQ} \\
  --modules-sub "${modules_sub}" \\
  --modules-all "${modules_all}" \\
${moduleSubnetworksArg}${renewalAssociationArg}${networkModesArg}  --layout-seed ${networkLayoutSeed} \\
  --layout-scale ${networkLayoutScale} \\
  --degree-scale ${networkDegreeScale} \\
  --degree-size-mode "${networkDegreeSizeMode}" \\
  --degree-min-area ${networkDegreeMinArea} \\
  --edge-width-scale ${networkEdgeWidthScale} \\
  --isa-scale ${networkIsaScale} \\
  --abundance-size-mode "${networkAbundanceSizeMode}" \\
  --abundance-reference ${networkAbundanceReference} \\
  --abundance-reference-area ${networkAbundanceReferenceArea} \\
  --abundance-min-area ${networkAbundanceMinArea} \\
  --abundance-max-area ${networkAbundanceMaxArea} \\
  --abundance-scale-power ${networkAbundanceScalePower} \\
  --max-labels ${networkMaxLabels}

ln -sf "${spieceasiOutputDirAbs}/network_modules_best_stats_all.tsv" network_modules_best_stats_all.tsv
ln -sf "${spieceasiOutputDirAbs}/spieceasi_network_layout_all.tsv" spieceasi_network_layout_all.tsv
if [[ "${microbialStateInterpretationEnabled}" == "true" ]]; then
  test -s "${spieceasiOutputDirAbs}/network_ecological_module_renewal_phase_POS_ALL.pdf"
  test -s "${spieceasiOutputDirAbs}/network_ecological_module_renewal_post_renewal_vs_stagnation_POS_ALL.pdf"
  test -s "${spieceasiOutputDirAbs}/network_asv_module_renewal_overlay.tsv"
fi
touch network.done
"""
}
process ASV_MAG_NETWORK {
    cpus 1
    conda "${asvMagNetworkCondaEnvPath}"

    input:
    path(graph, stageAs: 'asv_mag_network_graph.graphml')
    path(node_features)
    path(ecological_modules)
    path(ecological_module_selection)
    path(taxonomy_table)
    path(asv_counts)
    path(cruise_metadata)
    path(sample_metadata)
    path(isa_tables)
    path(dep_asv_mag, stageAs: 'dep_asv_mag.done')

    output:
    path("asv_mag_network.done"), emit: done
    path("analysis_eligible_mappings.tsv"), emit: analysis_mappings
    path("functional_modules.tsv"), emit: functional_modules, optional: true

    when:
    asvMagNetworkEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def magAbundanceSeqkitArg = asvMagNetworkMagAbundanceSeqkit ? """  --mag-abundance-seqkit "${asvMagNetworkMagAbundanceSeqkit}" \\
  --mag-abundance-seqkit-file-col "${asvMagNetworkMagAbundanceSeqkitFileCol}" \\
  --mag-abundance-seqkit-count-col "${asvMagNetworkMagAbundanceSeqkitCountCol}" \\
""" : ''
    def magAbundanceArg = asvMagNetworkMagAbundance ? """  --mag-abundance "${asvMagNetworkMagAbundance}" \\
${magAbundanceSeqkitArg}""" : ''
    def magTranscriptSeqkitArg = asvMagNetworkMagTranscriptAbundanceSeqkit ? """  --mag-transcript-abundance-seqkit "${asvMagNetworkMagTranscriptAbundanceSeqkit}" \\
  --mag-transcript-abundance-seqkit-file-col "${asvMagNetworkMagTranscriptAbundanceSeqkitFileCol}" \\
  --mag-transcript-abundance-seqkit-count-col "${asvMagNetworkMagTranscriptAbundanceSeqkitCountCol}" \\
""" : ''
    def magTranscriptAbundanceArg = asvMagNetworkMagTranscriptAbundance ? """  --mag-transcript-abundance "${asvMagNetworkMagTranscriptAbundance}" \\
  --mag-transcript-abundance-format "${asvMagNetworkMagTranscriptAbundanceFormat}" \\
  --mag-transcript-abundance-genome-col "${asvMagNetworkMagTranscriptAbundanceGenomeCol}" \\
  --mag-transcript-abundance-sample-col "${asvMagNetworkMagTranscriptAbundanceSampleCol}" \\
  --mag-transcript-abundance-value-col "${asvMagNetworkMagTranscriptAbundanceValueCol}" \\
  --mag-transcript-abundance-normalization "${asvMagNetworkMagTranscriptAbundanceNormalization}" \\
${magTranscriptSeqkitArg}""" : ''
    def functionalArgs = asvMagNetworkFunctionalAnnotations ? asvMagNetworkFunctionalAnnotations.collect { """  --functional-annotation "${it}" \\\n""" }.join('') : ''
    def asvTaxonomyConfidenceArg = asvMagLinkAsvTaxonomyMinConfidence != null ?
        """  --asv-taxonomy-min-confidence ${asvMagLinkAsvTaxonomyMinConfidence} \\\n""" : ''
    def displayAllEcologicalModulesArg = !networkModuleBestOnly ?
        """  --display-all-ecological-modules \\\n""" : ''
    """
set -euo pipefail
mkdir -p "${asvMagNetworkOutputDirAbs}"

echo "asv_mag_network.py md5: ${asvMagNetworkScriptHash}"
"\${CONDA_PREFIX}/bin/python" "${asvMagNetworkScriptPath}" \\
  --graph "${graph}" \\
  --node-features "${node_features}" \\
  --ecological-modules "${ecological_modules}" \\
  --ecological-module-selection "${ecological_module_selection}" \\
  --anchor-top-n ${asvMagNetworkAnchorTopN} \\
${displayAllEcologicalModulesArg}  --asv-mag-pairing "${asvMagLinkOutputDirAbs}/tables/asv2mag_pairing.tsv" \\
  --taxonomy "${taxonomy_table}" \\
  --asv-counts "${asv_counts}" \\
  --cruise-metadata "${cruise_metadata}" \\
  --sample-metadata "${sample_metadata}" \\
  --grouping-diagnostic-sample-col "${asvMagNetworkGroupingSampleCol}" \\
  --grouping-diagnostic-cruise-col "${asvMagNetworkGroupingCruiseCol}" \\
  --grouping-diagnostic-groups "${asvMagNetworkGroupingGroups}" \\
  --grouping-diagnostic-cruise-level-groups "${asvMagNetworkGroupingCruiseGroups}" \\
  --grouping-diagnostic-permutations ${asvMagNetworkGroupingPermutations} \\
  --isa-dir "${indicspeciesOutputDirAbs}" \\
  --biochemical-groupings "${asvMagNetworkBiochemicalGroupings}" \\
  --group-palette-map-json '${networkGroupPaletteJson}' \\
  --group-order-map-json '${networkGroupOrderJson}' \\
  --isa-q-threshold ${asvMagNetworkIsaQThreshold} \\
  --genome-summary "${asvMagLinkOutputDirAbs}/tables/asv2mag_genome_summary.tsv" \\
  --reference-catalog "${asvMagLinkOutputDirAbs}/references/barrnap_16s_reference_catalog.tsv" \\
  --outdir "${asvMagNetworkOutputDirAbs}" \\
  --prefix "${asvMagNetworkPrefix}" \\
  --asv-taxonomy-source "${asvMagNetworkAsvTaxonomySource}" \\
${asvTaxonomyConfidenceArg}  --mag-taxonomy-source "${asvMagNetworkMagTaxonomySource}" \\
  --mag-id-mode "${asvMagNetworkMagIdMode}" \\
  --mag-abundance-format "${asvMagNetworkMagAbundanceFormat}" \\
  --mag-abundance-genome-col "${asvMagNetworkMagAbundanceGenomeCol}" \\
  --mag-abundance-sample-col "${asvMagNetworkMagAbundanceSampleCol}" \\
  --mag-abundance-value-col "${asvMagNetworkMagAbundanceValueCol}" \\
  --min-shared-samples ${asvMagNetworkMinSharedSamples} \\
  --mag-knn ${asvMagNetworkMagKnn} \\
  --abundance-transform "${asvMagNetworkAbundanceTransform}" \\
  --mag-abundance-normalization "${asvMagNetworkMagAbundanceNormalization}" \\
  --functional-module-min-fraction ${asvMagNetworkFunctionalModuleMinFraction} \\
  --ambiguity-target-modules "${asvMagNetworkAmbiguityTargetModules}" \\
  --ambiguity-formats "${asvMagNetworkAmbiguityFormats}" \\
  --font-family "${asvMagNetworkFontFamily}" \\
  --min-pident ${asvMagNetworkMinPident} \\
${magAbundanceArg}${magTranscriptAbundanceArg}${functionalArgs}  --min-qcov ${asvMagNetworkMinQcov}
ln -sf "${asvMagNetworkOutputDirAbs}/mapping/${asvMagNetworkPrefix}_analysis_eligible_mappings.tsv" analysis_eligible_mappings.tsv
if [[ -f "${asvMagNetworkOutputDirAbs}/functional/${asvMagNetworkPrefix}_functional_modules.tsv" ]]; then
  ln -sf "${asvMagNetworkOutputDirAbs}/functional/${asvMagNetworkPrefix}_functional_modules.tsv" functional_modules.tsv
fi
touch asv_mag_network.done
"""
}

process ASV_MAG_CURTAINS {
    cpus 1
    conda "${asvMagCurtainsCondaEnvPath}"

    input:
    path(asv_mag_network_done)
    path(microbial_state_interpretation_done)

    output:
    path("asv_mag_curtains.done"), emit: done

    when:
    asvMagCurtainsEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def renewalArgs = asvMagCurtainsRenewalEvents ? """  --renewal-events "${asvMagCurtainsRenewalEvents}" \
  --renewal-date-col "${asvMagCurtainsRenewalDateCol}" \
""" : ''
    """
set -euo pipefail
echo "asv_mag_curtains.py md5: ${asvMagCurtainsScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"
export MPLCONFIGDIR="\$PWD/.matplotlib"
mkdir -p "\$MPLCONFIGDIR"
rm -rf asv_mag_curtains_staged

"\${CONDA_PREFIX}/bin/python" "${asvMagCurtainsScriptPath}" \
  --asv-mag-network-dir "${asvMagNetworkOutputDirAbs}" \
  --microbial-state-dir "${microbialStateInterpretationOutputDirAbs}" \
  --outdir asv_mag_curtains_staged \
  --prefix "${asvMagNetworkPrefix}" \
  --hybrid-palette "${microbialStateInterpretationHybridPalette}" \
  --hybrid-order "${microbialStateInterpretationHybridOrder}" \
  --mc-palette "${microbialStateInterpretationMcPalette}" \
  --mc-order "${microbialStateInterpretationMcOrder}" \
  --sample-col "${microbialStateInterpretationSampleCol}" \
  --date-col "${microbialStateInterpretationDateCol}" \
  --depth-col "${microbialStateInterpretationDepthCol}" \
  --hybrid-col "${microbialStateInterpretationHybridCol}" \
  --maximum-depth ${asvMagCurtainsMaximumDepth} \
  --minimum-point-area ${asvMagCurtainsMinimumPointArea} \
  --maximum-point-area ${asvMagCurtainsMaximumPointArea} \
  --asv-sample-point-area ${asvMagCurtainsAsvSamplePointArea} \
  --contrast-pseudocount-fpm ${asvMagCurtainsContrastPseudocountFpm} \
  --contrast-fold-threshold ${asvMagCurtainsContrastFoldThreshold} \
${renewalArgs}  --formats "${asvMagCurtainsFormats}"

"\${CONDA_PREFIX}/bin/python" -c 'from pathlib import Path; import shutil; src=Path("asv_mag_curtains_staged"); dst=Path("${asvMagCurtainsOutputDirAbs}"); shutil.rmtree(dst, ignore_errors=True); dst.parent.mkdir(parents=True, exist_ok=True); shutil.copytree(src, dst)'
[[ -f "${asvMagCurtainsOutputDirAbs}/tables/asv_mag_species_recruitment_curtain_manifest.tsv" ]] || { echo "Missing ASV-MAG curtain manifest" >&2; exit 1; }
[[ -f "${asvMagCurtainsOutputDirAbs}/tables/asv_mag_species_recruitment_curtain_points.tsv" ]] || { echo "Missing ASV-MAG curtain point table" >&2; exit 1; }
[[ -f "${asvMagCurtainsOutputDirAbs}/tables/asv_mag_species_recruitment_pairwise_contrast_summary.tsv" ]] || { echo "Missing ASV-MAG curtain contrast summary" >&2; exit 1; }
[[ -f "${asvMagCurtainsOutputDirAbs}/tables/asv_mag_curtain_asv_sample_anchors.tsv" ]] || { echo "Missing ASV-MAG curtain ASV sample anchors" >&2; exit 1; }
[[ -f "${asvMagCurtainsOutputDirAbs}/tables/asv_mag_curtain_renewal_onsets.tsv" ]] || { echo "Missing ASV-MAG curtain renewal-onset audit" >&2; exit 1; }
touch asv_mag_curtains.done
"""
}

process GROUP_GUILD_FUNCTION {
    cpus 1
    conda "${groupGuildFunctionCondaEnvPath}"

    input:
    path(metadata_table)
    path(asv_counts)
    path(modules_all)
    path(node_features)
    path(isa_tables)
    path(mappings)
    path(functional_modules)
    path(pca_scores)
    path(pca_explained)
    path(hybrid_assignments)
    path(hybrid_centroids)
    path(taxonomy_table)
    val(group_guild_function_script_hash)

    output:
    path("group_guild_function.done"), emit: done

    when:
    groupGuildFunctionEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def strictFontArg = groupGuildFunctionStrictFont ? '--strict-font' : ''
    def seasonValidation = groupGuildFunctionGroupings.split(',').collect { it.trim() }.contains('Season') ? """
[[ -s "${groupGuildFunctionOutputDirAbs}/tables/ecological_module_season_association.tsv" ]] || { echo "Missing ecological-module season association table" >&2; exit 1; }
[[ -s "${groupGuildFunctionOutputDirAbs}/tables/ecological_module_season_abundance.tsv" ]] || { echo "Missing ecological-module seasonal abundance table" >&2; exit 1; }
[[ -s "${groupGuildFunctionOutputDirAbs}/plots/ecological_module_seasonal_abundance_heatmap.pdf" ]] || { echo "Missing ecological-module seasonal heatmap" >&2; exit 1; }
""" : ''
    def pcaOverlayArg = groupGuildFunctionPcaOverlayEnabled ? """  --pca-scores "${pca_scores}" \\
  --pca-explained "${pca_explained}" \\
  --hybrid-assignments "${hybrid_assignments}" \\
  --hybrid-centroids "${hybrid_centroids}" \\
""" : ''
    """
set -euo pipefail
mkdir -p "${groupGuildFunctionOutputDirAbs}"
echo "group_guild_function.py md5: ${group_guild_function_script_hash}"
"\${CONDA_PREFIX}/bin/python" "${groupGuildFunctionScriptPath}" \\
  --metadata "${metadata_table}" \\
  --asv-counts "${asv_counts}" \\
  --modules "${modules_all}" \\
  --node-features "${node_features}" \\
  --isa-dir "${indicspeciesOutputDirAbs}" \\
  --accepted-mappings "${mappings}" \\
  --functional-modules "${functional_modules}" \\
  --taxonomy "${taxonomy_table}" \\
  --outdir "${groupGuildFunctionOutputDirAbs}" \\
  --groupings "${groupGuildFunctionGroupings}" \\
  --primary-groupings "${groupGuildFunctionPrimary}" \\
  --target-modules "${groupGuildFunctionTargets}" \\
  --target-module-min-fraction ${asvMagNetworkFunctionalModuleMinFraction} \\
  --sample-col "${groupGuildFunctionSampleCol}" \\
  --cruise-col "${groupGuildFunctionCruiseCol}" \\
  --permutations ${groupGuildFunctionPermutations} \\
  --seed ${groupGuildFunctionSeed} \\
  --max-modules ${groupGuildFunctionMaxModules} \\
  --heatmap-top-asvs ${groupGuildFunctionHeatmapTopAsvs} \\
  --pca-label-min-mean-abundance ${groupGuildFunctionPcaLabelMinMeanAbundancePct / 100.0} \\
${pcaOverlayArg}  --formats "${groupGuildFunctionFormats}" \\
  --font-family "${groupGuildFunctionFontFamily}" ${strictFontArg}
${seasonValidation}touch group_guild_function.done
"""
}

process ECOLOGICAL_CONTEXT_ATLAS {
    cpus 1
    conda "${ecologicalContextAtlasCondaEnvPath}"

    input:
    path(metadata_table)
    path(microbial_compartments_done)
    path(modules_all)
    path(node_features)
    path(network_layout)
    path(network_graph)
    path(taxonomy_table)
    path(titan_done)
    path(measurement_done)
    path(module_measurement_done)
    path(group_guild_done)
    path(mag_mappings)
    path(functional_modules)
    path(asv_mag_network_done)

    output:
    path("ecological_context_atlas.done"), emit: done

    when:
    ecologicalContextAtlasEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def rnaDnaRatioArg = ecologicalContextAtlasRnaDnaLog2TpmRatio ? """  --rna-dna-log2-tpm-ratio "${ecologicalContextAtlasRnaDnaLog2TpmRatio}" \\\n""" : ''
    def rnaDnaSeArg = ecologicalContextAtlasRnaDnaLog2TpmSe ? """  --rna-dna-log2-tpm-se "${ecologicalContextAtlasRnaDnaLog2TpmSe}" \\\n""" : ''
    """
set -euo pipefail
echo "ecological_context_atlas.py md5: ${ecologicalContextAtlasScriptHash}"
echo "asv_mag_network.py plotting md5: ${asvMagNetworkScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"

"\${CONDA_PREFIX}/bin/python" "${ecologicalContextAtlasScriptPath}" \
  --relative-abundance "${microbialCompartmentOutputDirAbs}/audit/asv_relative_abundance.tsv" \
  --metadata "${metadata_table}" \
  --microbial-compartments "${microbialCompartmentOutputDirAbs}/tables/microbial_compartments.tsv" \
  --taxonomy "${taxonomy_table}" \
  --modules "${modules_all}" \
  --node-features "${node_features}" \
  --network-layout "${network_layout}" \
  --network-graph "${network_graph}" \
  --titan-results "${titanOutputDirAbs}/tables/titan_taxon_results.tsv" \
  --correlations "${measurementAssociationOutputDirAbs}/tables/asv_measurement_spearman_long.tsv" \
  --module-compartment-associations "${groupGuildFunctionOutputDirAbs}/tables/compartment_matched_ecological_module_association.tsv" \
  --module-measurement-associations "${measurementAssociationOutputDirAbs}/tables/module_measurement_association.tsv" \
  --module-measurement-member-support "${measurementAssociationOutputDirAbs}/tables/module_measurement_member_support.tsv" \
  --indicator-dir "${indicspeciesOutputDirAbs}" \
  --mag-mappings "${mag_mappings}" \
  --functional-modules "${functional_modules}" \
  --target-functional-modules "${asvMagNetworkOutputDirAbs}/functional/${asvMagNetworkPrefix}_metagenome_ambiguous_mag_target_modules.tsv" \
  --heterogeneous-graph "${asvMagNetworkOutputDirAbs}/network/${asvMagNetworkPrefix}_metagenome_heterogeneous.graphml" \
  --heterogeneous-prefix "${asvMagNetworkPrefix}_metagenome" \
  --multiomics-dir "${asvMagNetworkOutputDirAbs}/abundance" \
  --outdir "${ecologicalContextAtlasOutputDirAbs}" \
  --sample-col "${ecologicalContextAtlasSampleCol}" \
  --cruise-col "${ecologicalContextAtlasCruiseCol}" \
  --depth-col "${ecologicalContextAtlasDepthCol}" \
  --date-col "${ecologicalContextAtlasDateCol}" \
  --context-cols "${ecologicalContextAtlasContextCols}" \
  --network-variables "${ecologicalContextAtlasNetworkVariables}" \
  --linkage-q-threshold ${ecologicalContextAtlasLinkageQThreshold} \
  --selected-asvs "${ecologicalContextAtlasSelectedAsvs}" \
  --context-sheet-source "${ecologicalContextAtlasContextSheetSource}" \
${rnaDnaRatioArg}${rnaDnaSeArg}  --rna-dna-max-log2-se ${ecologicalContextAtlasRnaDnaMaxLog2Se} \
  --rna-dna-mag-id-mode "${ecologicalContextAtlasRnaDnaMagIdMode}" \
  --formats "${ecologicalContextAtlasFormats}"

[[ -f "${ecologicalContextAtlasOutputDirAbs}/tables/asv_module_compartment_measurement_crosswalk.tsv" ]] || { echo "Missing complete ecological-linkage crosswalk" >&2; exit 1; }
[[ -f "${ecologicalContextAtlasOutputDirAbs}/tables/asv_module_compartment_measurement_significant.tsv" ]] || { echo "Missing significant ecological-linkage table" >&2; exit 1; }
[[ -f "${ecologicalContextAtlasOutputDirAbs}/tables/asv_module_compartment_measurement_mag_linked.tsv" ]] || { echo "Missing MAG-linked ecological-linkage table" >&2; exit 1; }
[[ -f "${ecologicalContextAtlasOutputDirAbs}/tables/asv_ecological_module_mag_crosswalk.tsv" ]] || { echo "Missing nonredundant ASV-module-MAG inventory" >&2; exit 1; }

touch ecological_context_atlas.done
"""
}

process MODULE_MAG_ANCHORS {
    cpus 1
    conda "${networkCondaEnvPath}"

    input:
    path(modules_all)
    path(node_features)
    path(taxonomy_table)
    path(asv_counts)
    path(metadata_table)
    path(asv_mag_mappings)
    path(dep_asv_mag)
    path(dep_graph_network)

    output:
    path("module_asv_anchor_table.tsv"), emit: asv_anchor_table
    path("module_mag_anchor_summary.tsv"), emit: module_summary
    path("sample_module_scores.tsv"), emit: sample_module_scores
    path("sample_top_modules.tsv"), emit: sample_top_modules
    path("sample_module_score_matrix.tsv"), emit: sample_module_matrix
    path("sample_module_score_heatmap.png"), optional: true, emit: sample_module_heatmap_png
    path("sample_module_score_heatmap.pdf"), optional: true, emit: sample_module_heatmap_pdf
    path("sample_module_score_heatmap.svg"), optional: true, emit: sample_module_heatmap_svg
    path("module_mag_anchors.done"), emit: done

    when:
    networkEnabled && asvMagLinkEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
mkdir -p "${spieceasiOutputDirAbs}"

"\${CONDA_PREFIX}/bin/python" "${moduleMagAnchorsScriptPath}" \\
  --modules "${modules_all}" \\
  --node-features "${node_features}" \\
  --taxonomy "${taxonomy_table}" \\
  --asv-mag-pairing "${asv_mag_mappings}" \\
  --asv-counts "${asv_counts}" \\
  --metadata "${metadata_table}" \\
  --sample-col "${metadataPlotsSampleCol}" \\
  --sample-code-col "${clustermapsSampleCodeCol}" \\
  --top-n ${asvMagNetworkAnchorTopN} \\
  --best-stats "${spieceasiOutputDirAbs}/network_modules_best_stats_all.tsv" \\
  --outdir "${spieceasiOutputDirAbs}"

ln -sf "${spieceasiOutputDirAbs}/module_asv_anchor_table.tsv" module_asv_anchor_table.tsv
ln -sf "${spieceasiOutputDirAbs}/module_mag_anchor_summary.tsv" module_mag_anchor_summary.tsv
ln -sf "${spieceasiOutputDirAbs}/sample_module_scores.tsv" sample_module_scores.tsv
ln -sf "${spieceasiOutputDirAbs}/sample_top_modules.tsv" sample_top_modules.tsv
ln -sf "${spieceasiOutputDirAbs}/sample_module_score_matrix.tsv" sample_module_score_matrix.tsv
for ext in png pdf svg; do
  if [[ -f "${spieceasiOutputDirAbs}/sample_module_score_heatmap.\${ext}" ]]; then
    ln -sf "${spieceasiOutputDirAbs}/sample_module_score_heatmap.\${ext}" "sample_module_score_heatmap.\${ext}"
  fi
done
touch module_mag_anchors.done
"""
}
process MASTER_SUMMARY {
    cpus 1
    conda "${masterSummaryCondaEnvPath}"

    input:
    path(asv_meta)
    path(asv_counts)
    path(dep_network, stageAs: 'dep_network.done')
    path(dep_sankey, stageAs: 'dep_sankey.done')
    path(dep_asv_mag, stageAs: 'dep_asv_mag.done')
    path(dep_optional, stageAs: 'dep_optional.done')

    output:
    path("ASV_master_long.tsv"), optional: true, emit: master_long
    path("ASV_master_count_wide.tsv"), optional: true, emit: master_count
    path("ASV_master_source_manifest.tsv"), optional: true, emit: master_manifest
    path("ASV_master_column_mapping.tsv"), optional: true, emit: master_colmap
    path("ASV_master_column_collisions_original.tsv"), optional: true, emit: master_collisions
    path("master_summary.done"), emit: done

    when:
    masterSummaryEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def whitelistArg = masterSummaryWhitelistCsv ? """  --whitelist "${masterSummaryWhitelistCsv}" \\\n""" : ''
    """
set -euo pipefail
mkdir -p "${masterSummaryOutputDirAbs}"

"\${CONDA_PREFIX}/bin/python" "${masterSummaryScriptPath}" \\
  --data-dir "${outputDir}" \\
  --asv-meta "${asv_meta}" \\
  --asv-counts "${asv_counts}" \\
  --clustermaps-dir "${masterSummaryClustermapsDirAbs}" \\
  --indicspecies-dir "${masterSummaryIndicspeciesDirAbs}" \\
  --spieceasi-dir "${masterSummarySpieceasiDirAbs}" \\
  --asv-mag-dir "${masterSummaryAsvMagDirAbs}" \\
${whitelistArg}  --outdir "${masterSummaryOutputDirAbs}" \\
  --max-direct-cols ${masterSummaryMaxDirectCols}

link_if_exists() {
  local src="\$1"
  local dest="\$2"
  if [[ -f "\${src}" ]]; then
    ln -sf "\${src}" "\${dest}"
  fi
}

link_if_exists "${masterSummaryOutputDirAbs}/ASV_master_long.tsv" "ASV_master_long.tsv"
link_if_exists "${masterSummaryOutputDirAbs}/ASV_master_count_wide.tsv" "ASV_master_count_wide.tsv"
link_if_exists "${masterSummaryOutputDirAbs}/ASV_master_source_manifest.tsv" "ASV_master_source_manifest.tsv"
link_if_exists "${masterSummaryOutputDirAbs}/ASV_master_column_mapping.tsv" "ASV_master_column_mapping.tsv"
link_if_exists "${masterSummaryOutputDirAbs}/ASV_master_column_collisions_original.tsv" "ASV_master_column_collisions_original.tsv"

touch master_summary.done
"""
}
process ASV_MAG_LINK {
    cpus asvMagLinkThreads
    conda "${asvMagLinkCondaEnvPath}"

    input:
    path(filtered_fasta)
    path(taxonomy_table)

    output:
    path("asv_mag_link.done"), emit: done
    path("asv2mag_pairing.tsv"), emit: pairing

    when:
    asvMagLinkEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def masterTsvArg = asvMagLinkMasterTsv ? """  --master-tsv "${asvMagLinkMasterTsv}" \\\n""" : ''
    def genomeDirArg = asvMagLinkGenomeDir ? """  --genome-fasta-dir "${asvMagLinkGenomeDir}" \\\n""" : ''
    def genomeQcDirArg = asvMagLinkGenomeQcDir ? """  --genome-qc-dir "${asvMagLinkGenomeQcDir}" \\\n""" : ''
    def genomeQcDirsArg = asvMagLinkGenomeQcDirs ? asvMagLinkGenomeQcDirs.collect { """  --genome-qc-dir "${it}" \\\n""" }.join('') : ''
    def idTokenIndexesArg = asvMagLinkIdTokenIndexes ? asvMagLinkIdTokenIndexes.collect { """  --id-token-index ${it} \\\n""" }.join('') : ''
    def barrnapDirArg = asvMagLinkBarrnapDir ? """  --barrnap-dir "${asvMagLinkBarrnapDir}" \\\n""" : ''
    def minCompletenessArg = asvMagLinkMinCompleteness != null ? """  --min-completeness ${asvMagLinkMinCompleteness} \\\n""" : ''
    def maxContaminationArg = asvMagLinkMaxContamination != null ? """  --max-contamination ${asvMagLinkMaxContamination} \\\n""" : ''
    def guncAssessmentArg = asvMagLinkGuncAssessmentValue ?
        """  --gunc-assessment-value "${asvMagLinkGuncAssessmentValue}" \\\n""" : ''
    def asvTaxonomyConfidenceArg = asvMagLinkAsvTaxonomyMinConfidence != null ?
        """  --asv-taxonomy "${taxonomy_table}" \\
  --asv-taxonomy-min-confidence ${asvMagLinkAsvTaxonomyMinConfidence} \\\n""" : ''
    def requireSpeciesArg = asvMagLinkRequireSpeciesAssignment ?
        """  --require-species-assignment \\\n""" : ''
    def minRrnaMarkerArg = asvMagLinkMinRrnaMarkerCount != null ?
        """  --min-rrna-marker-count ${asvMagLinkMinRrnaMarkerCount} \\\n""" : ''
    def autoBarrnapArg = asvMagLinkAutoBarrnap ?
        """  --auto-barrnap \\\n""" : ''
    """
set -euo pipefail
mkdir -p "${asvMagLinkOutputDirAbs}"

"\${CONDA_PREFIX}/bin/python" "${asvMagLinkScriptPath}" \\
  --asv-fasta "${filtered_fasta}" \\
${masterTsvArg}${barrnapDirArg}${genomeDirArg}${genomeQcDirArg}${genomeQcDirsArg}${idTokenIndexesArg}${autoBarrnapArg}  --outdir "${asvMagLinkOutputDirAbs}" \\
  --threads ${asvMagLinkThreads} \\
  --min-pident ${asvMagLinkMinPident} \\
  --min-qcov ${asvMagLinkMinQcov} \\
${asvTaxonomyConfidenceArg}${minCompletenessArg}${maxContaminationArg}${guncAssessmentArg}${requireSpeciesArg}${minRrnaMarkerArg}  --top-n ${asvMagLinkTopN}

"\${CONDA_PREFIX}/bin/python" "${plotAsvMagLinkScriptPath}" \\
  --input-dir "${asvMagLinkOutputDirAbs}" \\
  --top-n ${asvMagLinkPlotTopN}

ln -sf "${asvMagLinkOutputDirAbs}/tables/asv2mag_pairing.tsv" asv2mag_pairing.tsv
touch asv_mag_link.done
"""
}
