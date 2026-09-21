process SANKEY {
    cpus 1
    conda "${sankeyCondaEnvPath}"

    input:
    path(fastq_stats)
    path(filtered_stats)
    path(asv_counts)
    path(asv_decon_counts)
    path(asv_micro_counts)
    path(sample_manifest)

    output:
    path("sankey.done"), emit: done

    when:
    sankeyEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def keepTypesArg = sankeyKeepTypes && !sankeyKeepTypes.isEmpty() ? "  --keep-types \"${sankeyKeepTypes.join(',')}\" \\\n" : ''
    def verticalOrderArg = sankeyVerticalOrder && !sankeyVerticalOrder.isEmpty() ? "  --vertical-order \"${sankeyVerticalOrder.join(',')}\" \\\n" : ''
    def rawOutputPrefix = "${sankeyOutputPrefix}_raw"
    def labeledFlag = sankeyMakeLabeled ? "  --make-labeled \\\n" : ''
    def unlabeledFlag = sankeyMakeUnlabeled ? "  --make-unlabeled \\\n" : ''
"""
set -euo pipefail
echo "sankey_builder.py md5: ${sankeyScriptHash}"

"\${CONDA_PREFIX}/bin/python" "${sankeyScriptPath}" \\
  --data-dir "${outputDir}" \\
  --sub-dir "${sankeySubDir}" \\
  --metadata "${sankeyMetadataPath}" \\
  --sample-manifest "${sample_manifest}" \\
  --samp-col "${sankeySampCol}" \\
  --group1-col "${sankeyGroupCol}" \\
  --color-col "${sankeyColorCol}" \\
${keepTypesArg}${verticalOrderArg}  --fastq-stats "\${PWD}/${fastq_stats}" \\
  --filtered-stats "\${PWD}/${filtered_stats}" \\
  --asv-raw "\${PWD}/${asv_counts}" \\
  --asv-decon "\${PWD}/${asv_decon_counts}" \\
  --asv-micro "\${PWD}/${asv_micro_counts}" \\
  --title "${sankeyTitle}" \\
  --arrangement "${sankeyArrangement}" \\
  --output-prefix "${sankeyOutputPrefix}" \\
${labeledFlag}${unlabeledFlag}  --verbose

"\${CONDA_PREFIX}/bin/python" "${sankeyScriptPath}" \\
  --data-dir "${outputDir}" \\
  --sub-dir "${sankeySubDir}" \\
  --metadata "${sankeyMetadataPath}" \\
  --sample-manifest "${sample_manifest}" \\
  --samp-col "${sankeySampCol}" \\
  --group1-col "${sankeyGroupCol}" \\
  --color-col "${sankeyColorCol}" \\
${verticalOrderArg}  --fastq-stats "\${PWD}/${fastq_stats}" \\
  --filtered-stats "\${PWD}/${filtered_stats}" \\
  --asv-raw "\${PWD}/${asv_counts}" \\
  --asv-decon "\${PWD}/${asv_decon_counts}" \\
  --asv-micro "\${PWD}/${asv_micro_counts}" \\
  --title "${sankeyTitle}" \\
  --arrangement "${sankeyArrangement}" \\
  --output-prefix "${rawOutputPrefix}" \\
${labeledFlag}${unlabeledFlag}  --all-samples \\
  --verbose

touch sankey.done
"""
}
process GENERAL_STATS {
    cpus pipelineThreads
    conda "${generalStatsCondaEnvPath}"
    publishDir dirMap.stats, mode: 'copy', pattern: '*'

    input:
    path(concat_fasta)

    output:
    path("fastq_stats.tsv"), emit: fastq_stats
    path("fastp_fastqs.tsv"), emit: fastp_stats
    path("filtered_fastas.tsv"), emit: filtered_stats
    path("concat_fastas.tsv"), emit: concat_stats

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def rawArgs = generalStatsRawArgs
    def fastpArgs = generalStatsFastpArgs
    def filteredArgs = generalStatsFilteredArgs
    """
set -euo pipefail

run_seqkit() {
  local outfile="\$1"
  shift
  if [[ "\$#" -eq 0 ]]; then
    : > "\${outfile}"
    return
  fi
  seqkit stat -a -T -j ${task.cpus} -o "\${outfile}" "\$@"
}

run_seqkit fastq_stats.tsv ${rawArgs}
run_seqkit fastp_fastqs.tsv ${fastpArgs}
run_seqkit filtered_fastas.tsv ${filteredArgs}
run_seqkit concat_fastas.tsv "${concat_fasta}"
"""
}
process PLOT_METADATA {
    cpus pipelineThreads
    conda "${plotMetadataCondaEnvPath}"

    input:
    path(fastq_stats)
    path(asv_micro)
    path(asv_mito)
    path(taxonomy_table)
    path(sample_manifest)

    output:
    path("metadata_updated_micro.tsv"), emit: metadata_micro
    path("ASV_meta_micro.tsv"), emit: asv_meta_micro
    path("ASV_final.micro.tsv"), emit: asv_final_micro
    path("metadata_updated_mito.tsv"), optional: true, emit: metadata_mito
    path("ASV_meta_mito.tsv"), optional: true, emit: asv_meta_mito
    path("ASV_final.mito.tsv"), optional: true, emit: asv_final_mito

    when:
    metadataPlotsEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def includeRankAppend = metadataIncludeRank && !metadataIncludeRank.isEmpty() ?
        metadataIncludeRank.collect { "cmd+=( --include-rank \"${it}\" )" }.join('\n') : ''
    def metadataBiochemTable = metadataPlotsBiochemAssignmentsPath ?: ''
    def metadataBiochemIncludeCsv = metadataPlotsBiochemIncludeCols && !metadataPlotsBiochemIncludeCols.isEmpty() ? metadataPlotsBiochemIncludeCols.join(',') : ''
    def metadataBiochemMetaJoinCsv = metadataPlotsBiochemMetaJoinCols && !metadataPlotsBiochemMetaJoinCols.isEmpty() ? metadataPlotsBiochemMetaJoinCols.join(',') : ''
    def metadataBiochemJoinCsv = metadataPlotsBiochemJoinCols && !metadataPlotsBiochemJoinCols.isEmpty() ? metadataPlotsBiochemJoinCols.join(',') : ''
    def metadataStratTable = metadataPlotsStratificationTimeseriesPath ?: ''
    def metadataStratIncludeCsv = metadataPlotsStratIncludeCols && !metadataPlotsStratIncludeCols.isEmpty() ? metadataPlotsStratIncludeCols.join(',') : ''
    def metadataCruiseGroupTable = metadataPlotsCruiseGroupPath ?: ''
    def metadataCruiseGroupIncludeCsv = metadataPlotsCruiseGroupIncludeCols.join(',')
    def metadataKeepTypesCsv = metadataKeepTypes && !metadataKeepTypes.isEmpty() ? metadataKeepTypes.join(',') : ''
    def metadataSubtractionGroupsCsv = metadataPlotsSubtractionGroups && !metadataPlotsSubtractionGroups.isEmpty() ? metadataPlotsSubtractionGroups.join(',') : ''
    def metadataGroupOrderCsv = metadataPlotsGroupOrder && !metadataPlotsGroupOrder.isEmpty() ? metadataPlotsGroupOrder.join(',') : ''
    def metadataNormalizationColsCsv = metadataGroupNormalizationCols.join(',')
    def metadataMicroFile = "${outputDir}/metadata/metadata_updated_micro.tsv"
    def metadataMitoFile = "${outputDir}/mito/metadata/metadata_updated_mito.tsv"
    def asvMetaMicroFile = "${outputDir}/metadata/ASV_meta_micro.tsv"
    def asvMetaMitoFile = "${outputDir}/mito/metadata/ASV_meta_mito.tsv"
    def asvTargetMicroFile = "${outputDir}/ASVs/${filterCountsOutputName}"
    def asvTargetMitoFile = "${outputDir}/mito/ASVs/${filterCountsOutputName.replace('.tsv','.mito.tsv')}"
    def asvFinalMicroFile = "${outputDir}/ASVs/ASV_final.micro.tsv"
    def asvFinalMitoFile = "${outputDir}/mito/ASVs/ASV_final.mito.tsv"
    def asvTaxTable = "${outputDir}/taxonomy/ASV_SILVA_tax.full-length.vsearch.tsv"
"""
set -euo pipefail
echo "plot_metadata.py md5: ${plotMetadataScriptHash}"

cmd=(
  "\${CONDA_PREFIX}/bin/python" "${plotMetadataScriptPath}"
  --data-dir "${outputDir}"
  --sub-dir "${metadataPlotsSubDir}"
  --metadata "${metadataPlotsMetadataPath}"
  --taxonomy "${asvTaxTable}"
  --asv-micro "${asv_micro}"
  --asv-mito "${asv_mito}"
  --sample-id-col "${metadataPlotsSampleCol}"
  --group1-col "${metadataPlotsTypeCol}"
  --color-col "${metadataPlotsColorCol}"
  --subtraction-group-col "${metadataPlotsSubtractionCol}"
  --sample-manifest "${sample_manifest}"
  --make-micro
  --make-mito
  --verbose
)
cmd+=( --subtraction-groups "${metadataSubtractionGroupsCsv}" )
${includeRankAppend}
if [[ -n "${metadataKeepTypesCsv}" ]]; then
  cmd+=( --keep-types "${metadataKeepTypesCsv}" )
fi
if [[ -n "${metadataGroupOrderCsv}" ]]; then
  cmd+=( --group-order "${metadataGroupOrderCsv}" )
fi

if [[ -n "${metadataBiochemTable}" ]]; then
  if [[ -f "${metadataBiochemTable}" ]]; then
    cmd+=( --biochem-assignments "${metadataBiochemTable}" )
    cmd+=( --biochem-sample-col "${metadataPlotsBiochemSampleCol}" )
    if [[ -n "${metadataBiochemIncludeCsv}" ]]; then
      cmd+=( --biochem-include-cols "${metadataBiochemIncludeCsv}" )
    fi
    if [[ -n "${metadataBiochemMetaJoinCsv}" || -n "${metadataBiochemJoinCsv}" ]]; then
      if [[ -n "${metadataBiochemMetaJoinCsv}" && -n "${metadataBiochemJoinCsv}" ]]; then
        cmd+=( --biochem-meta-join-cols "${metadataBiochemMetaJoinCsv}" )
        cmd+=( --biochem-join-cols "${metadataBiochemJoinCsv}" )
      else
        echo "[w] Incomplete biochem join config; both metadata and biochem join col lists are required. Falling back to sample-id join." >&2
      fi
    fi
  else
    echo "[w] metadata biochem assignments table not found; skipping merge: ${metadataBiochemTable}" >&2
  fi
fi

if [[ -n "${metadataStratTable}" ]]; then
  if [[ -f "${metadataStratTable}" ]]; then
    cmd+=( --stratification-timeseries "${metadataStratTable}" )
    cmd+=( --stratification-meta-join-col "${metadataPlotsStratMetaJoinCol}" )
    cmd+=( --stratification-join-col "${metadataPlotsStratJoinCol}" )
    if [[ -n "${metadataStratIncludeCsv}" ]]; then
      cmd+=( --stratification-include-cols "${metadataStratIncludeCsv}" )
    fi
  else
    echo "[w] metadata stratification table not found; skipping merge: ${metadataStratTable}" >&2
  fi
fi

if [[ -n "${metadataCruiseGroupTable}" ]]; then
  if [[ -f "${metadataCruiseGroupTable}" ]]; then
    cmd+=( --cruise-group-assignments "${metadataCruiseGroupTable}" )
    cmd+=( --cruise-group-meta-join-col "${metadataPlotsCruiseGroupMetaJoinCol}" )
    cmd+=( --cruise-group-join-col "${metadataPlotsCruiseGroupJoinCol}" )
    cmd+=( --cruise-group-include-cols "${metadataCruiseGroupIncludeCsv}" )
  else
    echo "[w] cruise group assignments table not found; skipping merge: ${metadataCruiseGroupTable}" >&2
  fi
fi

if [[ "${metadataGroupNormalizationEnabled}" == "true" ]]; then
  cmd+=( --normalize-group-cols "${metadataNormalizationColsCsv}" )
  cmd+=( --normalize-group-pattern '${metadataGroupNormalizationPattern}' )
  cmd+=( --normalize-group-replacement "${metadataGroupNormalizationReplacement}" )
  if [[ "${metadataGroupNormalizationPreserveSource}" == "true" ]]; then
    cmd+=( --preserve-normalized-source )
  fi
fi

"\${cmd[@]}"

link_if_exists() {
  local src="\$1"
  local dest="\$2"
  if [[ -f "\${src}" ]]; then
    ln -sf "\${src}" "\${dest}"
  fi
}

link_if_exists "${metadataMicroFile}" "metadata_updated_micro.tsv"
link_if_exists "${asvMetaMicroFile}" "ASV_meta_micro.tsv"
link_if_exists "${asvFinalMicroFile}" "ASV_final.micro.tsv"
link_if_exists "${metadataMitoFile}" "metadata_updated_mito.tsv"
link_if_exists "${asvMetaMitoFile}" "ASV_meta_mito.tsv"
link_if_exists "${asvFinalMitoFile}" "ASV_final.mito.tsv"
"""
}

process ASV_TIME_DEPTH_CURTAIN {
    cpus 1
    conda "${asvTimeDepthCurtainCondaEnvPath}"

    input:
    path(metadata_table)

    output:
    path("asv_time_depth_curtain.done"), emit: done

    when:
    asvTimeDepthCurtainEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def curtainFormatChecks = asvTimeDepthCurtainFormats.split(',').collect { fmt ->
        "[[ -f \"${asvTimeDepthCurtainOutputDirAbs}/plots/asv_hybrid_compartment_time_depth_curtain.${fmt.trim()}\" ]] || { echo \"Missing ASV curtain ${fmt.trim()}\" >&2; exit 1; }"
    }.join('\n')
    def curtainOptionalArgs = []
    if( asvTimeDepthCurtainSourceMode == 'basin' ) {
        curtainOptionalArgs << "--basin-grid \"${asvTimeDepthCurtainBasinGrid}\""
        curtainOptionalArgs << "--basin-cells \"${asvTimeDepthCurtainBasinCells}\""
    }
    if( asvTimeDepthCurtainRenewalEvents ) {
        curtainOptionalArgs << "--renewal-events \"${asvTimeDepthCurtainRenewalEvents}\""
        curtainOptionalArgs << "--renewal-date-col \"${asvTimeDepthCurtainRenewalDateCol}\""
    }
    def curtainOptionalArgText = curtainOptionalArgs ?
        " \\\n  " + curtainOptionalArgs.join(" \\\n  ") : ''
    """
set -euo pipefail
echo "asv_time_depth_curtain.py md5: ${asvTimeDepthCurtainScriptHash}"
echo "aspire cache generation: ${aspireCacheGeneration}"

"\${CONDA_PREFIX}/bin/python" "${asvTimeDepthCurtainScriptPath}" \\
  --metadata "${metadata_table}" \\
  --outdir "${asvTimeDepthCurtainOutputDirAbs}" \\
  --sample-col "${asvTimeDepthCurtainSampleCol}" \\
  --cruise-col "${asvTimeDepthCurtainCruiseCol}" \\
  --date-col "${asvTimeDepthCurtainDateCol}" \\
  --depth-col "${asvTimeDepthCurtainDepthCol}" \\
  --group-col "${asvTimeDepthCurtainGroupCol}" \\
  --palette '${asvTimeDepthCurtainPalette}' \\
  --group-order '${asvTimeDepthCurtainOrder}' \\
  --exclude-label-pattern '${asvTimeDepthCurtainExcludePattern}' \\
  --maximum-depth ${asvTimeDepthCurtainMaximumDepth} \\
  --depth-step ${asvTimeDepthCurtainDepthStep} \\
  --time-subdivisions-per-month ${asvTimeDepthCurtainTimeSubdivisions} \\
  --time-sigma-months ${asvTimeDepthCurtainTimeSigma} \\
  --contour-visual-depth-sigma-m ${asvTimeDepthCurtainContourDepthSigma} \\
  --base-point-size ${asvTimeDepthCurtainBasePointSize} \\
  --asv-point-size ${asvTimeDepthCurtainAsvPointSize} \\
  --asv-point-edge-width ${asvTimeDepthCurtainAsvPointEdgeWidth} \\
  --formats "${asvTimeDepthCurtainFormats}"${curtainOptionalArgText}

[[ -f "${asvTimeDepthCurtainOutputDirAbs}/tables/asv_hybrid_compartment_time_depth_grid.tsv" ]] || { echo "Missing ASV curtain grid" >&2; exit 1; }
[[ -f "${asvTimeDepthCurtainOutputDirAbs}/tables/asv_hybrid_compartment_time_depth_renewals.tsv" ]] || { echo "Missing ASV curtain renewal-onset audit" >&2; exit 1; }
if [[ "${asvTimeDepthCurtainSourceMode}" == "basin" ]]; then
  [[ -f "${asvTimeDepthCurtainOutputDirAbs}/tables/asv_basin_time_depth_point_matches.tsv" ]] || { echo "Missing BASIN-ASV curtain point audit" >&2; exit 1; }
fi
${curtainFormatChecks}
touch asv_time_depth_curtain.done
"""
}

process PLOT_UPSET {
    cpus pipelineThreads
    conda "${plotUpsetCondaEnvPath}"

    input:
    path(metadata_table)

    output:
    path("plot_upset.done"), emit: done

    when:
    plotUpsetEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def taxonomyArg = plotUpsetTaxonomyPath ? """  --taxonomy-path "${plotUpsetTaxonomyPath}" \\\n""" : ''
    def groupOrderArg = plotUpsetGroupOrder && !plotUpsetGroupOrder.isEmpty() ? """  --group-order "${plotUpsetGroupOrder.join(',')}" \\\n""" : ''
    def subsetGroupsArg = plotUpsetSubsetGroups && !plotUpsetSubsetGroups.isEmpty() ? """  --subset-groups "${plotUpsetSubsetGroups.join(',')}" \\\n""" : ''
    def vennSubsetGroupsArg = plotUpsetVennSubsetGroups && !plotUpsetVennSubsetGroups.isEmpty() ? """  --venn-subset-groups "${plotUpsetVennSubsetGroups.join(',')}" \\\n""" : ''
    def groupPaletteArg = plotUpsetGroupPalette ? """  --group-palette "${plotUpsetGroupPalette}" \\\n""" : ''
    def skipVennArg = plotUpsetSkipVenn ? "  --skip-venn \\\n" : ''
    def rawOnlyArg = plotUpsetRawOnly ? "  --raw-only \\\n" : ''
    def finalOnlyArg = plotUpsetFinalOnly ? "  --final-only \\\n" : ''
    def rawMicroMetadataPath = "${outputDir}/metadata/metadata_updated_micro_raw.tsv"
    def rawMicroFinalAsvPath = "${outputDir}/ASVs/ASV_final_raw.micro.tsv"
    def rawMicroAsvTargetPath = "${outputDir}/ASVs/ASV_target.micro.tsv"
    def rawMitoMetadataPath = "${outputDir}/mito/metadata/metadata_updated_mito_raw.tsv"
    def rawMitoFinalAsvPath = "${outputDir}/mito/ASVs/ASV_final_raw.mito.tsv"
    def rawMitoAsvTargetPath = "${outputDir}/mito/ASVs/ASV_target.mito.tsv"
    def rawMetadataPathSingle = plotUpsetDomain == 'mito' ? rawMitoMetadataPath : rawMicroMetadataPath
    def rawFinalAsvPathSingle = plotUpsetDomain == 'mito' ? rawMitoFinalAsvPath : rawMicroFinalAsvPath
    def rawAsvTargetPathSingle = plotUpsetDomain == 'mito' ? rawMitoAsvTargetPath : rawMicroAsvTargetPath
    """
set -euo pipefail

"\${CONDA_PREFIX}/bin/python" "${plotUpsetScriptPath}" \\
  --data-dir "${outputDir}" \\
  --subdir "${plotUpsetSubDir}" \\
  --domain "${plotUpsetDomain}" \\
${taxonomyArg}  --sample-id-col "${plotUpsetSampleIdCol}" \\
  --group-col "${plotUpsetGroupCol}" \\
  --color-col "${plotUpsetColorCol}" \\
${groupPaletteArg}${groupOrderArg}${subsetGroupsArg}${vennSubsetGroupsArg}${skipVennArg}${rawOnlyArg}${finalOnlyArg}  --formats "${plotUpsetFormats}" \\
  --font-size ${plotUpsetFontSize} \\
  --max-intersections ${plotUpsetMaxIntersections}

if [[ "${plotUpsetDomain}" == "both" ]]; then
  "\${CONDA_PREFIX}/bin/python" "${plotUpsetScriptPath}" \\
    --data-dir "${outputDir}" \\
    --subdir "${plotUpsetSubDir}" \\
    --domain "micro" \\
${taxonomyArg}    --sample-id-col "${plotUpsetSampleIdCol}" \\
    --group-col "${plotUpsetGroupCol}" \\
    --color-col "${plotUpsetColorCol}" \\
${groupPaletteArg}${groupOrderArg}${subsetGroupsArg}${vennSubsetGroupsArg}${skipVennArg}${rawOnlyArg}${finalOnlyArg}    --formats "${plotUpsetFormats}" \\
    --font-size ${plotUpsetFontSize} \\
    --max-intersections ${plotUpsetMaxIntersections} \\
    --metadata-path "${rawMicroMetadataPath}" \\
    --asv-raw-path "${rawMicroAsvTargetPath}" \\
    --asv-final-path "${rawMicroFinalAsvPath}" \\
    --output-tag raw

  "\${CONDA_PREFIX}/bin/python" "${plotUpsetScriptPath}" \\
    --data-dir "${outputDir}" \\
    --subdir "${plotUpsetSubDir}" \\
    --domain "mito" \\
${taxonomyArg}    --sample-id-col "${plotUpsetSampleIdCol}" \\
    --group-col "${plotUpsetGroupCol}" \\
    --color-col "${plotUpsetColorCol}" \\
${groupPaletteArg}${groupOrderArg}${subsetGroupsArg}${vennSubsetGroupsArg}${skipVennArg}${rawOnlyArg}${finalOnlyArg}    --formats "${plotUpsetFormats}" \\
    --font-size ${plotUpsetFontSize} \\
    --max-intersections ${plotUpsetMaxIntersections} \\
    --metadata-path "${rawMitoMetadataPath}" \\
    --asv-raw-path "${rawMitoAsvTargetPath}" \\
    --asv-final-path "${rawMitoFinalAsvPath}" \\
    --output-tag raw
else
  "\${CONDA_PREFIX}/bin/python" "${plotUpsetScriptPath}" \\
    --data-dir "${outputDir}" \\
    --subdir "${plotUpsetSubDir}" \\
    --domain "${plotUpsetDomain}" \\
${taxonomyArg}    --sample-id-col "${plotUpsetSampleIdCol}" \\
    --group-col "${plotUpsetGroupCol}" \\
    --color-col "${plotUpsetColorCol}" \\
${groupPaletteArg}${groupOrderArg}${subsetGroupsArg}${vennSubsetGroupsArg}${skipVennArg}${rawOnlyArg}${finalOnlyArg}    --formats "${plotUpsetFormats}" \\
    --font-size ${plotUpsetFontSize} \\
    --max-intersections ${plotUpsetMaxIntersections} \\
    --metadata-path "${rawMetadataPathSingle}" \\
    --asv-raw-path "${rawAsvTargetPathSingle}" \\
    --asv-final-path "${rawFinalAsvPathSingle}" \\
    --output-tag raw
fi

touch plot_upset.done
"""
}
process BUBBLEPLOTTER {
    cpus pipelineThreads
    conda "${bubbleplotterCondaEnvPath}"

    input:
    path(asv_meta)

    output:
    path("bubbleplotter.done"), emit: done

    when:
    bubbleplotterEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def noAutoSizeArg = bubbleplotterNoAutoSize ? "  --no-auto-size \\\n" : ''
    def bubbleplotterGroup1OrderArg = bubbleplotterGroup1Order && !bubbleplotterGroup1Order.isEmpty() ? """  --group1-order "${bubbleplotterGroup1Order.join(',')}" \\\n""" : ''
    def bubbleplotterGroup2OrderArg = bubbleplotterGroup2Order && !bubbleplotterGroup2Order.isEmpty() ? """  --group2-order "${bubbleplotterGroup2Order.join(',')}" \\\n""" : ''
    """
set -euo pipefail
mkdir -p "${bubbleplotterOutputDirAbs}"

"\${CONDA_PREFIX}/bin/python" "${bubbleplotterScriptPath}" \\
  --input "${asv_meta}" \\
  --output-prefix "${bubbleplotterOutputPrefixAbs}" \\
  --count-col "${bubbleplotterCountCol}" \\
  --sample-col "${bubbleplotterSampleCol}" \\
  --group1-col "${bubbleplotterDepthCol}" \\
  --color-col "${bubbleplotterColorCol}" \\
  --group2-col "${bubbleplotterMonthCol}" \\
${bubbleplotterGroup1OrderArg}${bubbleplotterGroup2OrderArg}${noAutoSizeArg}  --formats "${bubbleplotterFormats}" \\
  --figsize "${bubbleplotterFigsize}" \\
  --bubble-scale ${bubbleplotterScale}

touch bubbleplotter.done
"""
}
process UMAP_CLUSTERING {
    cpus pipelineThreads
    conda "${umapClusteringCondaEnvPath}"

    input:
    path(asv_meta)

    output:
    path("umap_clustering.done"), emit: done

    when:
    umapClusteringEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def umapGroup1OrderArg = umapClusteringGroup1Order && !umapClusteringGroup1Order.isEmpty() ? """  --group1-order "${umapClusteringGroup1Order.join(',')}" \\\n""" : ''
    def umapGroup2OrderArg = umapClusteringGroup2Order && !umapClusteringGroup2Order.isEmpty() ? """  --group2-order "${umapClusteringGroup2Order.join(',')}" \\\n""" : ''
    def umapNoScaleArg = umapClusteringNoScale ? "  --no-scale \\\n" : ''
    """
set -euo pipefail
mkdir -p "${umapClusteringOutputDirAbs}"

"\${CONDA_PREFIX}/bin/python" "${umapClusteringScriptPath}" \\
  --input "${asv_meta}" \\
  --output-prefix "${umapClusteringOutputPrefixAbs}" \\
  --count-col "${umapClusteringCountCol}" \\
  --sample-col "${umapClusteringSampleCol}" \\
  --group1-col "${umapClusteringDepthCol}" \\
  --color-col "${umapClusteringColorCol}" \\
  --group2-col "${umapClusteringSecondaryCol}" \\
  --group1-palette "${umapClusteringGroup1Palette}" \\
  --group2-palette "${umapClusteringGroup2Palette}" \\
${umapGroup1OrderArg}${umapGroup2OrderArg}  --formats "${umapClusteringFormats}" \\
  --normalize "${umapClusteringNormalize}" \\
  --transform "${umapClusteringTransform}" \\
  --n-neighbors ${umapClusteringNeighbors} \\
  --min-dist ${umapClusteringMinDist} \\
  --umap-metric "${umapClusteringMetric}" \\
  --min-cluster-size ${umapClusteringMinClusterSize} \\
  --min-samples ${umapClusteringMinSamples} \\
  --hdbscan-metric "${umapClusteringHdbscanMetric}" \\
${umapNoScaleArg}  --random-state ${umapClusteringRandomState}

touch umap_clustering.done
"""
}
process ASV_BATCH_CORRECTION {
    cpus pipelineThreads
    conda "${batchCorrectionCondaEnvPath}"

    input:
    path(metadata_table)
    path(asv_meta)
    path(asv_counts)

    output:
    path("asv_clr_selected.tsv"), emit: asv_clr_selected
    path("asv_clr_after_correction.tsv"), emit: asv_clr_after
    path("asv_corrected_abundance.features_rows.tsv"), emit: asv_corrected_counts
    path("asv_corrected_pseudocount.features_rows.tsv"), emit: asv_corrected_counts_int
    path("asv_selected_abundance.features_rows.tsv"), emit: asv_selected_counts
    path("asv_selected_pseudocount.features_rows.tsv"), emit: asv_selected_counts_int
    path("batch_correction_decision.tsv"), emit: correction_decision
    path("batch_correction_countspace_preservation.png"), emit: countspace_plot
    path("batch_correction_countspace_preservation_metrics.tsv"), emit: countspace_metrics
    path("batch_correction_umap_comparison.png"), emit: umap_plot
    path("batch_correction_statistics.tsv"), emit: correction_stats
    path("umap_hdbscan_results.tsv"), emit: umap_results

    when:
    batchCorrectionEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def bioCovArg = batchBiologicalCovariates ? """  --biological-covariates "${batchBiologicalCovariates}" \\\n""" : ''
    def minSamplesArg = batchHdbscanMinSamples != null ? "  --hdbscan-min-samples ${batchHdbscanMinSamples} \\\n" : ''
    def optimizeFlag = batchOptimize ? "  --optimize-clustering \\\n" : ''
    def conqurBatchRefArg = batchConqurBatchRef ? """  --conqur-batch-ref "${batchConqurBatchRef}" \\\n""" : ''
    def conqurLogisticLassoFlag = batchConqurLogisticLasso ? "  --conqur-logistic-lasso \\\n" : ''
    def conqurSimpleMatchFlag = batchConqurSimpleMatch ? "  --conqur-simple-match \\\n" : ''
    def conqurInterpltFlag = batchConqurInterplt ? "  --conqur-interplt \\\n" : ''
    def conqurAutoInstallFlag = batchConqurAutoInstall ? "  --conqur-auto-install \\\n" : ''
    def asvClrAfterFile = "${batchCorrectionOutputDirAbs}/asv_clr_after_correction.tsv"
    def asvCorrectedFeaturesFile = "${batchCorrectionOutputDirAbs}/asv_corrected_abundance.features_rows.tsv"
    def asvCorrectedPseudoFeaturesFile = "${batchCorrectionOutputDirAbs}/asv_corrected_pseudocount.features_rows.tsv"
    def asvSelectedClrFile = "${batchCorrectionOutputDirAbs}/asv_clr_selected.tsv"
    def asvSelectedFeaturesFile = "${batchCorrectionOutputDirAbs}/asv_selected_abundance.features_rows.tsv"
    def asvSelectedPseudoFeaturesFile = "${batchCorrectionOutputDirAbs}/asv_selected_pseudocount.features_rows.tsv"
    def batchCorrectionDecisionFile = "${batchCorrectionOutputDirAbs}/batch_correction_decision.tsv"
    def countspacePlotFile = "${batchCorrectionOutputDirAbs}/batch_correction_countspace_preservation.png"
    def countspaceMetricsFile = "${batchCorrectionOutputDirAbs}/batch_correction_countspace_preservation_metrics.tsv"
    def umapComparisonPngFile = "${batchCorrectionOutputDirAbs}/batch_correction_umap_comparison.png"
    def batchCorrectionStatsFile = "${batchCorrectionOutputDirAbs}/batch_correction_statistics.tsv"
    def umapResultsFile = "${batchCorrectionOutputDirAbs}/umap_hdbscan_results.tsv"
    """
set -euo pipefail

"\${CONDA_PREFIX}/bin/python" "${batchCorrectionScriptPath}" \\
  --data-dir "${outputDir}" \\
  --asv "\$PWD/${asv_counts}" \\
  --metadata "\$PWD/${metadata_table}" \\
  --asv-meta "\$PWD/${asv_meta}" \\
  --sample-id-col "${batchCorrectionSampleIdCol}" \\
  --batch-col "${batchCorrectionBatchCol}" \\
  --output-dir "${batchCorrectionOutputDir}" \\
  --asv-orientation "${batchCorrectionOrientation}" \\
  --conqur-mode "${batchConqurMode}" \\
  --correction-policy "${batchCorrectionPolicy}" \\
  --auto-min-sample-rho ${batchAutoMinSampleRho} \\
  --auto-min-bray-rho ${batchAutoMinBrayRho} \\
  --auto-max-batch-eta-ratio ${batchAutoMaxBatchEtaRatio} \\
  --auto-min-batch-eta-drop ${batchAutoMinBatchEtaDrop} \\
  --auto-min-bio-eta-ratio ${batchAutoMinBioEtaRatio} \\
  --conqur-num-core ${batchConqurNumCore} \\
${conqurBatchRefArg}${conqurLogisticLassoFlag}  --conqur-quantile-type "${batchConqurQuantileType}" \\
${conqurSimpleMatchFlag}  --conqur-lambda-quantile "${batchConqurLambdaQuantile}" \\
${conqurInterpltFlag}  --conqur-delta ${batchConqurDelta} \\
${conqurAutoInstallFlag}${bioCovArg}  --umap-neighbors ${batchUmapNeighbors} \\
  --umap-min-dist ${batchUmapMinDist} \\
  --hdbscan-min-cluster-size ${batchHdbscanMinClusterSize} \\
${minSamplesArg}  --hdbscan-selection-method "${batchHdbscanSelectionMethod}" \\
  --target-clusters "${batchTargetClusters}" \\
  --n-features-plot ${batchNFeaturesPlot} \\
  --biological-color-col "${batchBiologicalColorCols}" \\
  --color-palette-col "${batchColorPaletteCols}" \\
  --biological-palettes-json '${batchBiologicalPalettesJson}' \\
  --random-state ${batchRandomState} \\
${optimizeFlag}  --verbose

if [[ ! -f "${asvClrAfterFile}" ]]; then
  echo "Missing batch correction output: ${asvClrAfterFile}" >&2
  exit 1
fi
ln -sf "${asvClrAfterFile}" asv_clr_after_correction.tsv
if [[ ! -f "${asvCorrectedFeaturesFile}" ]]; then
  echo "Missing batch correction output: ${asvCorrectedFeaturesFile}" >&2
  exit 1
fi
ln -sf "${asvCorrectedFeaturesFile}" asv_corrected_abundance.features_rows.tsv
if [[ ! -f "${asvCorrectedPseudoFeaturesFile}" ]]; then
  echo "Missing batch correction output: ${asvCorrectedPseudoFeaturesFile}" >&2
  exit 1
fi
ln -sf "${asvCorrectedPseudoFeaturesFile}" asv_corrected_pseudocount.features_rows.tsv
if [[ ! -f "${asvSelectedClrFile}" ]]; then
  echo "Missing batch correction selected output: ${asvSelectedClrFile}" >&2
  exit 1
fi
ln -sf "${asvSelectedClrFile}" asv_clr_selected.tsv
if [[ ! -f "${asvSelectedFeaturesFile}" ]]; then
  echo "Missing batch correction selected output: ${asvSelectedFeaturesFile}" >&2
  exit 1
fi
ln -sf "${asvSelectedFeaturesFile}" asv_selected_abundance.features_rows.tsv
if [[ ! -f "${asvSelectedPseudoFeaturesFile}" ]]; then
  echo "Missing batch correction selected output: ${asvSelectedPseudoFeaturesFile}" >&2
  exit 1
fi
ln -sf "${asvSelectedPseudoFeaturesFile}" asv_selected_pseudocount.features_rows.tsv
if [[ ! -f "${batchCorrectionDecisionFile}" ]]; then
  echo "Missing batch correction decision output: ${batchCorrectionDecisionFile}" >&2
  exit 1
fi
ln -sf "${batchCorrectionDecisionFile}" batch_correction_decision.tsv
if [[ ! -f "${countspacePlotFile}" ]]; then
  echo "Missing batch correction output: ${countspacePlotFile}" >&2
  exit 1
fi
ln -sf "${countspacePlotFile}" batch_correction_countspace_preservation.png
if [[ ! -f "${countspaceMetricsFile}" ]]; then
  echo "Missing batch correction output: ${countspaceMetricsFile}" >&2
  exit 1
fi
ln -sf "${countspaceMetricsFile}" batch_correction_countspace_preservation_metrics.tsv
if [[ ! -f "${umapComparisonPngFile}" ]]; then
  echo "Missing batch correction output: ${umapComparisonPngFile}" >&2
  exit 1
fi
ln -sf "${umapComparisonPngFile}" batch_correction_umap_comparison.png
if [[ ! -f "${batchCorrectionStatsFile}" ]]; then
  echo "Missing batch correction output: ${batchCorrectionStatsFile}" >&2
  exit 1
fi
ln -sf "${batchCorrectionStatsFile}" batch_correction_statistics.tsv
if [[ ! -f "${umapResultsFile}" ]]; then
  echo "Missing batch correction output: ${umapResultsFile}" >&2
  exit 1
fi
ln -sf "${umapResultsFile}" umap_hdbscan_results.tsv
"""
}
process ASV_META_FROM_CORRECTED {
    cpus 1
    conda "${batchCorrectionCondaEnvPath}"

    input:
    path(asv_meta)
    path(corrected_counts)

    output:
    path("ASV_meta_micro.corrected.tsv"), emit: asv_meta_corrected

    when:
    batchCorrectionEnabled && (
        bubbleplotterEnabled || umapClusteringEnabled || clustermapsEnabled ||
        vocCorrelationEnabled || measurementAssociationEnabled || masterSummaryEnabled ||
        powerAnalysisEnabled || taxonomyPatientAwareEnabled || lungStatusAnalysisEnabled
    )

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail

"\${CONDA_PREFIX}/bin/python" - "${asv_meta}" "${corrected_counts}" "ASV_meta_micro.corrected.tsv" <<'PY'
import sys
import pandas as pd

asv_meta_path, corrected_counts_path, out_path = sys.argv[1:4]
sample_col = "${batchCorrectionSampleIdCol}"
asv_col = "ASV_ID"
count_col = "count"
count_alias_col = "corr_count"

meta = pd.read_csv(asv_meta_path, sep='\\t')
if sample_col not in meta.columns or asv_col not in meta.columns:
    raise ValueError(f"Expected columns '{sample_col}' and '{asv_col}' in {asv_meta_path}")

corr = pd.read_csv(corrected_counts_path, sep='\\t', index_col=0)
corr.index = corr.index.astype(str)
corr.columns = corr.columns.astype(str)

corr_long = corr.stack().rename(count_col).reset_index()
corr_long.columns = [asv_col, sample_col, count_col]
corr_long[count_col] = pd.to_numeric(corr_long[count_col], errors='coerce').fillna(0.0).clip(lower=0.0)

candidate_cols = [c for c in meta.columns if c not in {sample_col, asv_col, count_col, count_alias_col}]
asv_only_cols = []
sample_only_cols = []
for col in candidate_cols:
    asv_n = meta.groupby(asv_col, dropna=False)[col].nunique(dropna=False).max()
    sample_n = meta.groupby(sample_col, dropna=False)[col].nunique(dropna=False).max()
    if asv_n <= 1 and sample_n > 1:
        asv_only_cols.append(col)
    else:
        sample_only_cols.append(col)

sample_meta = meta[[sample_col] + sample_only_cols].drop_duplicates(subset=[sample_col])
asv_meta_df = meta[[asv_col] + asv_only_cols].drop_duplicates(subset=[asv_col])

out = corr_long.merge(sample_meta, on=sample_col, how='left')
out = out.merge(asv_meta_df, on=asv_col, how='left')
out = out[out[count_col] > 0]
out[count_alias_col] = out[count_col]

front = [sample_col, asv_col, count_col, count_alias_col]
remaining = [c for c in out.columns if c not in front]
out = out[front + remaining]
out.to_csv(out_path, sep='\\t', index=False)

print(f"[i] Wrote corrected ASV meta table: {out_path}")
print(f"[i] Rows: {len(out)}")
PY
"""
}
process OUTLIER_CHECKER {
    cpus pipelineThreads
    conda "${outlierCondaEnvPath}"

    input:
    path(asv_clr)
    path(metadata_table)

    when:
    outlierEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def groupColsArg = outlierGroupCols.join(',')
    def isoFlag = outlierUseIso ? "  --use-iso \\\n" : ''
    def svmFlag = outlierUseSvm ? "  --use-svm \\\n" : ''
    def hdbFlag = outlierUseHdb ? "  --use-hdb \\\n" : ''
    def preTransFlag = outlierPreTransformed ? "  --pre-transformed \\\n" : ''
    def scaleFlag = outlierScale ? "  --scale \\\n" : ''
    def hdbMinSamplesArg = outlierHdbMinSamples != null ? "  --hdbscan-min-samples ${outlierHdbMinSamples} \\\n" : ''
    def asvClrAfterFile = "${batchCorrectionOutputDirAbs}/asv_clr_after_correction.tsv"


    """
set -euo pipefail

"\${CONDA_PREFIX}/bin/python" "${outlierCheckerScriptPath}" \\
  --data-dir "${outputDir}" \\
  --asv "${asvClrAfterFile}" \\
  --metadata "\$PWD/${metadata_table}" \\
  --sample-id-col "${outlierSampleIdCol}" \\
  --output-dir "${outlierOutputDirAbs}" \\
  --group-cols "${groupColsArg}" \\
  --asv-orientation "${outlierOrientation}" \\
  --transform "${outlierTransform}" \\
${preTransFlag}${scaleFlag}${isoFlag}${svmFlag}${hdbFlag}  --vote-threshold ${outlierVoteThreshold} \\
  --iso-contamination "${outlierIsoContamination}" \\
  --iso-estimators ${outlierIsoEstimators} \\
  --iso-random-state ${outlierIsoRandomState} \\
  --svm-kernel "${outlierSvmKernel}" \\
  --svm-gamma "${outlierSvmGamma}" \\
  --svm-nu ${outlierSvmNu} \\
  --hdbscan-min-cluster-size ${outlierHdbMinClusterSize} \\
${hdbMinSamplesArg}  --hdbscan-metric "${outlierHdbMetric}" \\
  --verbose
"""
}
process COLLECTORS_CURVE {
    cpus pipelineThreads
    conda "${collectorsCondaEnvPath}"

    input:
    path(asv_counts)
    path(metadata_table)

    when:
    collectorsEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def collectorsGroupOrderArg = collectorsGroupOrder && !collectorsGroupOrder.isEmpty() ? """  --group-order "${collectorsGroupOrder.join(',')}" \\\n""" : ''
    def collectorsGroupColorsArg = collectorsGroupColors ? """  --group-colors "${collectorsGroupColors}" \\\n""" : ''
    """
set -euo pipefail

"\${CONDA_PREFIX}/bin/python" "${collectorsCurveScriptPath}" \\
  --counts "${asv_counts}" \\
  --meta "${metadata_table}" \\
  --sample-col "${collectorsSampleCol}" \\
  --group-col "${collectorsGroupCol}" \\
  --color-col "${collectorsColorCol}" \\
${collectorsGroupColorsArg}${collectorsGroupOrderArg}  --permutations ${collectorsPermutations} \\
  --seed ${collectorsSeed} \\
  --out_prefix "${collectorsOutPrefixAbs}" \\
  --title "${collectorsTitle}" \\
  --formats "${collectorsFormats}" \\
  --xpad ${collectorsXpad} \\
  --max-cols ${collectorsMaxCols} \\
  --show-perms ${collectorsShowPerms} \\
  --presence-threshold ${collectorsPresenceThreshold}
"""
}
process DIVERSITY_ANALYSIS {
    cpus pipelineThreads
    conda "${diversityCondaEnvPath}"

    input:
    path(metadata_table)
    path(asv_counts)
    path(matched_cohort, stageAs: "matched_cohort/*")

    output:
    path("diversity.done"), emit: done

    when:
    diversityEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def secondaryColArg = diversitySecondaryCol ? """  --secondary-col "${diversitySecondaryCol}" \\\n""" : ''
    def excludeGroupsArg = diversityExcludeGroups && !diversityExcludeGroups.isEmpty() ? """  --exclude-groups "${diversityExcludeGroups.join(',')}" \\\n""" : ''
    def groupOrderArg = diversityGroupOrder && !diversityGroupOrder.isEmpty() ? """  --group-order "${diversityGroupOrder.join(',')}" \\\n""" : ''
    def diversityBlockArg = diversityBlockCol ? """  --block-col "${diversityBlockCol}" \\\n""" : ''
    def verboseFlag = diversityVerbose ? "  --verbose\n" : ''
    def diversityRunMitoFlag = diversityRunMito ? '1' : '0'
    def diversityPatientAwareEnabledFlag = diversityPatientAwareEnabled ? '1' : '0'
    def diversityPatientAwareExcludeContralateralFlag = diversityPatientAwareExcludeContralateral ? 'TRUE' : 'FALSE'
    def diversityPatientAwareRequireCompleteTypesFlag = diversityPatientAwareRequireCompleteTypes ? '1' : '0'
    def cohortArg = matched_cohort ? """  --cohort-table "${matched_cohort}" \\\n""" : ''
    """
set -euo pipefail
mkdir -p "${diversityOutputDirAbs}"
mkdir -p "${diversityMitoOutputDirAbs}"

if [[ "${diversityRunMitoFlag}" == "1" && -f "${diversityMitoInputPath}" ]]; then
  "\${CONDA_PREFIX}/bin/python" "${calcDivScriptPath}" \\
    --micro-table "${asv_counts}" \\
    --mito-table "${diversityMitoInputPath}" \\
    --outdir "${diversityOutputDirAbs}" \\
    --mito-outdir "${diversityMitoOutputDirAbs}"
else
  "\${CONDA_PREFIX}/bin/python" "${calcDivScriptPath}" \\
    --micro-table "${asv_counts}" \\
    --outdir "${diversityOutputDirAbs}"
fi

"\${CONDA_PREFIX}/bin/python" "${plotDiversityScriptPath}" \\
  --metadata "${metadata_table}" \\
  --sample-col "${diversitySampleCol}" \\
  --group-col "${diversityGroupCol}" \\
  --color-col "${diversityColorCol}" \\
  --group-palette "${diversityGroupPalette}" \\
  --secondary-palette "${diversitySecondaryPalette}" \\
${cohortArg}${secondaryColArg}${excludeGroupsArg}${groupOrderArg}  --alpha-table "${diversityOutputDirAbs}/shannon.tsv" \\
  --distance-bray "${diversityOutputDirAbs}/bray.tsv" \\
  --distance-jaccard "${diversityOutputDirAbs}/jaccard.tsv" \\
  --output-dir "${diversityOutputDirAbs}" \\
  --umap-neighbors ${diversityUmapNeighbors} \\
  --umap-min-dist ${diversityUmapMinDist} \\
${diversityBlockArg}  --permanova-perms ${diversityPermutations} \\
  --random-state ${diversityRandomState} \\
${verboseFlag}

if [[ "${diversityRunMitoFlag}" == "1" && -f "${diversityMitoOutputDirAbs}/shannon.mito.tsv" && -f "${diversityMitoOutputDirAbs}/bray.mito.tsv" && -f "${diversityMitoOutputDirAbs}/jaccard.mito.tsv" ]]; then
  "\${CONDA_PREFIX}/bin/python" "${plotDiversityScriptPath}" \\
    --metadata "${metadata_table}" \\
    --sample-col "${diversitySampleCol}" \\
    --group-col "${diversityGroupCol}" \\
    --color-col "${diversityColorCol}" \\
    --group-palette "${diversityGroupPalette}" \\
    --secondary-palette "${diversitySecondaryPalette}" \\
${cohortArg}${secondaryColArg}${excludeGroupsArg}${groupOrderArg}    --alpha-table "${diversityOutputDirAbs}/shannon.tsv" \\
    --distance-bray "${diversityOutputDirAbs}/bray.tsv" \\
    --distance-jaccard "${diversityOutputDirAbs}/jaccard.tsv" \\
    --output-dir "${diversityOutputDirAbs}" \\
    --mito-mode \\
    --mito-alpha "${diversityMitoOutputDirAbs}/shannon.mito.tsv" \\
    --mito-bray "${diversityMitoOutputDirAbs}/bray.mito.tsv" \\
    --mito-jaccard "${diversityMitoOutputDirAbs}/jaccard.mito.tsv" \\
    --mito-output-dir "${diversityMitoOutputDirAbs}" \\
    --umap-neighbors ${diversityUmapNeighbors} \\
    --umap-min-dist ${diversityUmapMinDist} \\
${diversityBlockArg}    --permanova-perms ${diversityPermutations} \\
    --random-state ${diversityRandomState} \\
${verboseFlag}
fi

if [[ "${diversityPatientAwareEnabledFlag}" == "1" ]]; then
  mkdir -p "${diversityPatientAwareOutputDirAbs}"

  DIVERSITY_PATIENT_AWARE_ARGS=()
  if [[ "${diversityPatientAwareRequireCompleteTypesFlag}" == "1" ]]; then
    DIVERSITY_PATIENT_AWARE_ARGS+=(--require-complete-types)
  fi

  Rscript "${brayPatientAwareScriptPath}" \\
    --data-wide "${asv_counts}" \\
    --data-long "${metadata_table}" \\
    --sample-col "${diversityPatientAwareSampleCol}" \\
    --patient-col "${diversityPatientAwarePatientCol}" \\
    --case-col "${diversityPatientAwareCaseCol}" \\
    --type-col "${diversityPatientAwareTypeCol}" \\
    --sample-types "${diversityPatientAwareSampleTypes}" \\
    --contralateral-col "${diversityPatientAwareContralateralCol}" \\
    --cancer-site-col "${diversityPatientAwareCancerSiteCol}" \\
    --lung-side-col "${diversityPatientAwareLungSideCol}" \\
    --contralateral-value "${diversityPatientAwareContralateralValue}" \\
    --contralateral-sample-types "${diversityPatientAwareContralateralTypes}" \\
    --exclude-contralateral-in-cancer ${diversityPatientAwareExcludeContralateralFlag} \\
    --transform "${diversityPatientAwareTransform}" \\
    --permutations ${diversityPatientAwarePermutations} \\
    --seed ${diversityPatientAwareSeed} \\
    --outdir "${diversityPatientAwareOutputDirAbs}" \\
    "\${DIVERSITY_PATIENT_AWARE_ARGS[@]}"

  "\${CONDA_PREFIX}/bin/python" "${plotBrayPatientAwareScriptPath}" \\
    --indir "${diversityPatientAwareOutputDirAbs}" \\
    --outdir "${diversityPatientAwareOutputDirAbs}/figures"
fi

touch diversity.done
"""
}
