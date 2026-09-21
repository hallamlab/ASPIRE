process FASTP_QC {
    tag { meta.sample_id }
    cpus sampleThreads
    maxForks maxParallelSampleTasks
    conda "${condaEnvPath}"
    publishDir dirMap.fastp, mode: 'copy', pattern: '*', saveAs: { filename ->
        if( filename == 'R1.fastq.gz' ) {
            return "${meta.sample_id}_R1.fastq.gz"
        }
        if( filename == 'R2.fastq.gz' ) {
            return "${meta.sample_id}_R2.fastq.gz"
        }
        if( filename == 'fastp.json' ) {
            return "${meta.sample_id}.fastp.json"
        }
        if( filename == 'fastp.html' ) {
            return "${meta.sample_id}.fastp.html"
        }
        return filename
    }

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("R1.fastq.gz"), path("R2.fastq.gz"), emit: reads
    path("fastp.json"), emit: fastp_json
    path("fastp.html"), emit: fastp_html

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def trimFrontR1 = fastpTrimValues.front_r1
    def trimTailR1  = fastpTrimValues.tail_r1
    def trimFrontR2 = fastpTrimValues.front_r2
    def trimTailR2  = fastpTrimValues.tail_r2
    if( meta.paired && r2 ) {
        return """
fastp \\
  -i "${r1}" -I "${r2}" \\
  -o R1.fastq.gz \\
  -O R2.fastq.gz \\
  -f ${trimFrontR1} -t ${trimTailR1} \\
  -F ${trimFrontR2} -T ${trimTailR2} \\
  -j fastp.json \\
  -h fastp.html \\
  -w ${task.cpus}
"""
    }
    return """
fastp \\
  -i "${r1}" \\
  -o R1.fastq.gz \\
  -f ${trimFrontR1} -t ${trimTailR1} \\
  -j fastp.json \\
  -h fastp.html \\
  -w ${task.cpus}
\nln -sf R1.fastq.gz R2.fastq.gz\n
"""
}
process MERGE_READS {
    tag { meta.sample_id }
    cpus sampleThreads
    maxForks maxParallelSampleTasks
    conda "${condaEnvPath}"
    publishDir dirMap.merge, mode: 'copy', saveAs: { filename ->
        filename == 'merged.fastq.gz' ? "${meta.sample_id}.merged.fastq.gz" : filename
    }

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("merged.fastq.gz")

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def allowStagger = mergeAllowStagger ? '--fastq_allowmergestagger' : ''
    if( meta.paired && r2 ) {
        """
set -euo pipefail
vsearch --fastq_mergepairs "${r1}" \\
        --reverse "${r2}" \\
        --fastqout merged.fastq \\
        --fastq_maxdiffs ${mergeMaxDiffs} \\
        --fastq_minovlen ${mergeMinOverlap} \\
        --fastq_truncqual ${mergeTruncQuality} \\
        ${allowStagger} \\
        --threads ${task.cpus}
gzip -n merged.fastq
"""
    } else {
        """
set -euo pipefail
if [[ "${r1}" == *.gz ]]; then
  gunzip -c "${r1}" > merged.fastq
else
  cat "${r1}" > merged.fastq
fi
gzip -n merged.fastq
"""
    }
}
process FILTER_READS {
    tag { meta.sample_id }
    cpus sampleThreads
    maxForks maxParallelSampleTasks
    conda "${condaEnvPath}"
    publishDir dirMap.filter, mode: 'copy', saveAs: { filename ->
        filename == 'filtered.fasta.gz' ? "${meta.sample_id}.filtered.fasta.gz" : filename
    }

    input:
    tuple val(meta), path(merged_fastq)

    output:
    tuple val(meta), path("filtered.fasta.gz")

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
gzip -cd "${merged_fastq}" > merged.fastq
vsearch --fastx_filter merged.fastq \\
        --fastq_maxee ${filterMaxEe} \\
        --fastq_minlen ${filterMinLen} \\
        --fastq_maxlen ${filterMaxLen} \\
        --fastaout filtered.fasta
gzip -n filtered.fasta
"""
}
process RELABEL_FILTERED {
    tag { meta.sample_id }
    cpus 1
    maxForks maxParallelSampleTasks
    conda "${condaEnvPath}"
    publishDir dirMap.concat, mode: 'copy', pattern: '*.fasta.gz'

    input:
    tuple val(meta), path(filtered_fasta)

    output:
    tuple val(meta), path("*.filtered.relabel.fasta.gz"), emit: relabeled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def labelSep = concatLabelSep
    def relabeledOut = "${meta.sample_id}.filtered.relabel.fasta.gz"
    """
awk -v pref="${meta.sample_id}" -v sep="${labelSep}" '{
  if (\$0 ~ /^>/) {
    sub(/^>[^:]*:/, ">" pref sep, \$0)
    if (\$0 !~ "^>" pref sep) \$0 = ">" pref sep substr(\$0, 2)
  }
  print
}' <(gzip -cd "${filtered_fasta}") | gzip -n > "${relabeledOut}"
"""
}
process CONCAT_FASTAS {
    cpus 1
    conda "${condaEnvPath}"
    publishDir dirMap.concat, mode: 'copy', pattern: '*.fasta.gz'

    input:
    path(filtered_fastas)

    output:
    path("concat.fasta.gz"), emit: concat_for_derep
    path("concat_counts.fasta.gz"), emit: concat_for_counts

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def concatInputs = filtered_fastas.collect { "\"${it}\"" }.join(' ')
    """
set -euo pipefail
for f in ${concatInputs}; do
  gzip -cd "\${f}" || cat "\${f}"
done > concat.fasta
gzip -n -c concat.fasta > concat.fasta.gz
cp concat.fasta.gz concat_counts.fasta.gz
"""
}
process DEREPLICATE {
    cpus pipelineThreads
    conda "${condaEnvPath}"
    publishDir dirMap.derep, mode: 'copy', pattern: '*'

    input:
    path(concat_fasta)

    output:
    path("derep.fasta.gz")

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
gzip -cd "${concat_fasta}" > concat.fasta
vsearch --derep_fulllength concat.fasta \\
        --output derep.fasta \\
        --sizeout \\
        --threads ${task.cpus} \\
        --log "${dirMap.logs}/derep.log"
gzip -n derep.fasta
"""
}
process SINA_TRIM {
    cpus sinaThreads
    conda "${sinaCondaEnvPath}"
    publishDir dirMap.sina, mode: 'copy', pattern: '*'

    input:
    path(derep_fasta)

    output:
    path("derep_trimmed.fasta.gz"), emit: trimmed_fasta
    path("derep_SINA.fasta.gz")
    path("derep_SINA.log")
    path("derep_v_regions.tsv")

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def parseVerboseArg = sinaVerbose ? ' --verbose' : ''
    def trimTargetArg = sinaTrimTarget ? " -t \"${sinaTrimTarget}\"" : ''
    def keepGapsArg = sinaKeepGaps ? ' --keep-gaps' : ''
    """
set -euo pipefail
gzip -cd "${derep_fasta}" > derep_input.fasta
sina \\
    -i derep_input.fasta \\
    -o derep_SINA.fasta \\
    -r "${sinaReferencePath}" \\
    -v \\
    -p ${task.cpus} \\
    --log-file derep_SINA.log

"\${CONDA_PREFIX}/bin/python" "${parseSinaScriptPath}" \\
  --log derep_SINA.log \\
  --output derep_v_regions.tsv${parseVerboseArg}

"\${CONDA_PREFIX}/bin/python" "${trimSinaScriptPath}" \\
  -m derep_v_regions.tsv \\
  -f derep_SINA.fasta \\
  -r "${sinaRegionsArg}"${trimTargetArg} \\
  -o derep_trimmed.fasta \\
  --id-column "${sinaIdColumn}" \\
  --threads ${task.cpus} \\
  --batch-size ${sinaBatchSize}${keepGapsArg}
gzip -n derep_SINA.fasta
gzip -n derep_trimmed.fasta
"""
}
process DENOISE {
    cpus pipelineThreads
    conda "${condaEnvPath}"
    publishDir dirMap.denoise, mode: 'copy', pattern: '*'

    input:
    path(trimmed_fasta)

    output:
    path("centroids.fasta.gz")

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def unoiseCfg = config.unoise ?: [:]
    """
set -euo pipefail
gzip -cd "${trimmed_fasta}" > derep_trimmed.fasta
vsearch --cluster_unoise derep_trimmed.fasta \\
        --centroids centroids.fasta \\
        --sizein --sizeout --relabel ASV \\
        --minsize ${unoiseCfg.min_size != null ? unoiseCfg.min_size : 8} \\
        --threads ${task.cpus} \\
        --log "${dirMap.logs}/denoise.log"
gzip -n centroids.fasta
"""
}
process CHIMERA_CHECK {
    cpus pipelineThreads
    conda "${condaEnvPath}"
    publishDir dirMap.nochi, mode: 'copy', pattern: '*'

    input:
    path(centroids)

    output:
    path("nochimeras.fasta.gz")

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
gzip -cd "${centroids}" > centroids.fasta
vsearch --uchime3_denovo centroids.fasta \\
        --nonchimeras nochimeras.fasta \\
        --sizein \\
        --threads ${task.cpus} \\
        --log "${dirMap.logs}/nochimera.log"
gzip -n nochimeras.fasta
"""
}
process CREATE_COUNT_MATRIX {
    cpus pipelineThreads
    conda "${condaEnvPath}"
    publishDir dirMap.asv, mode: 'copy', pattern: '*'

    input:
    path(concat_fasta)
    path(nochimeras)

    output:
    tuple path("ASV_counts.tsv"), path("ASVs.fasta.gz"), emit: count_matrix

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def mappingCfg = config.read_mapping ?: [:]
    def mappingIdentity = mappingCfg.min_identity != null ? (mappingCfg.min_identity as double) : 0.999d
    if( mappingIdentity <= 0 || mappingIdentity > 1 ) {
        error 'read_mapping.min_identity must be greater than 0 and at most 1'
    }
    """
set -euo pipefail
gzip -cd "${nochimeras}" > ASVs.fasta
gzip -cd "${concat_fasta}" > concat_counts.fasta
vsearch --usearch_global concat_counts.fasta \\
        --db ASVs.fasta \\
        --id ${mappingIdentity} \\
        --otutabout ASV_counts.tsv \\
        --threads ${task.cpus} \\
        --log "${dirMap.logs}/count.log"
gzip -n ASVs.fasta
"""
}
process FILTER_TABLE {
    conda "${condaEnvPath}"
    publishDir dirMap.asv, mode: 'copy', pattern: '*'

    input:
    tuple path(count_table), path(asv_fasta)

    output:
    tuple path("ASV_filtered.tsv"), path("ASVs_filtered.fasta.gz"), emit: filtered

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def tableCfg = config.table_filter ?: [:]
    """
set -euo pipefail
gzip -cd "${asv_fasta}" > ASVs.fasta
"\${CONDA_PREFIX}/bin/python" "${tableScriptFile}" \\
       "${count_table}" \\
       ASV_filtered.tsv \\
       ${tableCfg.min_sample_sum != null ? tableCfg.min_sample_sum : 5000} \\
       ${tableCfg.min_asv_sum != null ? tableCfg.min_asv_sum : 0.01} \\
       ASVs.fasta \\
       ASVs_filtered.fasta
gzip -n ASVs_filtered.fasta
"""
}

process PREPROCESS_DATASET {
    tag "canonical preprocessing dataset"
    cpus 1
    conda "${condaEnvPath}"
    publishDir publicOutputDir, mode: 'copy', pattern: 'preprocessing_dataset', overwrite: true

    input:
    path(asv_counts)
    path(asv_filtered)
    path(asv_sequences)
    path(asv_target)
    path(asv_target_micro)
    path(asv_target_mito)
    path(asv_target_decon)
    path(taxonomy)
    path(fastq_stats)
    path(fastp_stats)
    path(filtered_stats)
    path(concat_stats)
    path(sample_manifest)
    val(preprocessing_parameters)

    output:
    path("preprocessing_dataset"), emit: dataset

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
"\${CONDA_PREFIX}/bin/python" "${preprocessDatasetBuildScriptPath}" \\
  --outdir preprocessing_dataset \\
  --preprocessing-parameters '${preprocessing_parameters.replace("'", "'\\''")}' \\
  --asv-counts "${asv_counts}" \\
  --asv-filtered "${asv_filtered}" \\
  --asv-sequences "${asv_sequences}" \\
  --asv-target "${asv_target}" \\
  --asv-target-micro "${asv_target_micro}" \\
  --asv-target-mito "${asv_target_mito}" \\
  --asv-target-decon "${asv_target_decon}" \\
  --taxonomy "${taxonomy}" \\
  --fastq-stats "${fastq_stats}" \\
  --fastp-stats "${fastp_stats}" \\
  --filtered-stats "${filtered_stats}" \\
  --concat-stats "${concat_stats}" \\
  --sample-manifest "${sample_manifest}"
"""
}

process VALIDATE_PREPROCESS_DATASET {
    tag "canonical preprocessing dataset"
    cpus 1
    // The all-phase workflow receives this directory from a Nextflow work
    // path, whereas analysis-only receives the published canonical copy.
    // Deep hashing makes those byte-identical directories share one cache
    // identity and prevents a validation-path change from cascading through
    // every downstream analysis task.
    cache 'deep'
    conda "${condaEnvPath}"

    input:
    path(dataset)

    output:
    path("validated_preprocessing_dataset"), emit: dataset

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
"\${CONDA_PREFIX}/bin/python" "${preprocessDatasetValidateScriptPath}" "${dataset}"
cp -aL "${dataset}" validated_preprocessing_dataset
"""
}
