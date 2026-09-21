process TAXONOMY {
    cpus taxonomyThreads
    conda "${taxonomyCondaEnvPath}"
    publishDir dirMap.taxonomy, mode: 'copy', pattern: '*'

    input:
    path(filtered_fasta)

    output:
    path("${taxonomyUppercaseGzName}")
    path("${taxonomyOutputName}"), emit: taxonomy_table
    path("${taxonomyStatsName}")

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    """
set -euo pipefail
gzip -cd "${filtered_fasta}" > filtered_input.fasta
"\${CONDA_PREFIX}/bin/python" - <<'PY' filtered_input.fasta '${taxonomyUppercasePlainName}'
import sys
from pathlib import Path
src = Path(sys.argv[1])
dst = Path(sys.argv[2])
with src.open() as inp, dst.open('w') as out:
    for line in inp:
        if line.startswith('>'):
            out.write(line)
        else:
            out.write(line.strip().upper() + '\\n')
PY

"\${CONDA_PREFIX}/bin/python" "${taxonomyScriptPath}" \\
  --input-fasta "${taxonomyUppercasePlainName}" \\
  --ref-taxonomy "${taxonomyRefTaxonomy}" \\
  --ref-seqs "${taxonomyRefSequences}" \\
  --output-tsv "${taxonomyOutputName}" \\
  --stats-output "${taxonomyStatsName}" \\
  --threads ${task.cpus}
gzip -n -c "${taxonomyUppercasePlainName}" > "${taxonomyUppercaseGzName}"
"""
}
process PREPARE_BLAST_DATABASES {
    cpus 1
    conda "${mitomasterCondaEnvPath}"
    publishDir dirMap.reference, mode: 'copy', pattern: 'blast_databases'

    input:
    val(mito_source)
    val(contaminant_source)

    output:
    tuple path("blast_databases/mitochondrial"), path("blast_databases/contaminants"), emit: databases

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def mitoSourceCommand = mitoBlastFastaPath ?
        (mitoBlastFastaPath.toLowerCase().endsWith('.gz') ?
            "gzip -cd \"${mitoBlastFastaPath}\" > blast_databases/mitochondrial/source.fasta" :
            "cp \"${mitoBlastFastaPath}\" blast_databases/mitochondrial/source.fasta") :
        "blastdbcmd -db \"${mitoBlastDbPath}\" -entry all -out blast_databases/mitochondrial/source.fasta"
    def contaminantSourceCommand = mitoBiofFastaPath ?
        (mitoBiofFastaPath.toLowerCase().endsWith('.gz') ?
            "gzip -cd \"${mitoBiofFastaPath}\" > blast_databases/contaminants/source.fasta" :
            "cp \"${mitoBiofFastaPath}\" blast_databases/contaminants/source.fasta") :
        "blastdbcmd -db \"${mitoBiofDbPath}\" -entry all -out blast_databases/contaminants/source.fasta"
    """
set -euo pipefail
mkdir -p blast_databases/mitochondrial blast_databases/contaminants
${mitoSourceCommand}
${contaminantSourceCommand}
makeblastdb -in blast_databases/mitochondrial/source.fasta -dbtype nucl -parse_seqids -out blast_databases/mitochondrial/db
makeblastdb -in blast_databases/contaminants/source.fasta -dbtype nucl -parse_seqids -out blast_databases/contaminants/db
"""
}
process MITOMASTER {
    cpus mitoBlastThreads
    conda "${mitomasterCondaEnvPath}"
    publishDir mitoOutputDirPath, mode: 'copy', pattern: '*'

    input:
    tuple path(filtered_table), path(filtered_fasta)
    tuple path(mito_db_dir), path(contaminant_db_dir)

    output:
    tuple path("mitomaster_output.tsv"), path("mito_ncbi.blast6.tsv"), path("ssu_pipeline_contaminants.blast6.tsv"), emit: mito_artifacts

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def filteredFastaPath = filtered_fasta.toString().trim()
    def mitomasterCommand = mitoRunMitomaster ? """
"\${CONDA_PREFIX}/bin/python" "${mitomasterScriptPath}" \\
  --data-dir "${mitoChunkDirPath}" \\
  --glob-pattern "*.fa*" \\
  --recursive \\
  --output-file mitomaster_output.tsv \\
  --max-workers ${mitomasterWorkers} \\
  --timeout ${mitomasterTimeout} \\
  --retries ${mitomasterRetries} \\
  --header-mode "${mitomasterHeaderMode}" \\
  --overwrite
""" : "printf 'Sequence_ID\\thaplo\\n' > mitomaster_output.tsv"
    """
set -euo pipefail
rm -rf "${mitoChunkDirPath}"
mkdir -p "${mitoChunkDirPath}"
gzip -cd "${filteredFastaPath}" > filtered_input.fasta
FILTERED_FASTA=\$(realpath filtered_input.fasta)

seqkit split -s ${mitoChunkSize} -O "${mitoChunkDirPath}" "\${FILTERED_FASTA}"

${mitomasterCommand}

blastn -query "\${FILTERED_FASTA}" \\
  -db "${mito_db_dir}/db" \\
  -outfmt "6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore" \\
  -out mito_ncbi.blast6.tsv \\
  -num_threads ${task.cpus}

blastn -query "\${FILTERED_FASTA}" \\
  -db "${contaminant_db_dir}/db" \\
  -outfmt "6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore" \\
  -out ssu_pipeline_contaminants.blast6.tsv \\
  -num_threads ${task.cpus}
"""
}
process MITO_DECONTAM {
    cpus mitoBlastThreads
    conda "${mitoCheckerCondaEnvPath}"
    publishDir mitoOutputDirPath, mode: 'copy', pattern: '*'

    input:
    tuple path(mitomaster_file), path(mito_blast), path(biof_blast)
    path(taxonomy_table)

    output:
    path("${mitoPrefix}.master.tsv"), emit: nontarget_table
    path("${mitoPrefix}.summary_*.tsv"), optional: true, emit: summary_tables
    path("${mitoPrefix}_*.svg"), optional: true, emit: svg_plots
    path("${mitoPrefix}_*.pdf"), optional: true, emit: pdf_plots
    path("${mitoPrefix}_*.png"), optional: true, emit: png_plots

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def noPlotsFlag = mitoNoPlots ? ' --no-plots' : ''
    """
"\${CONDA_PREFIX}/bin/python" "${mitoCheckerScriptPath}" \\
  --mitomaster-file "${mitomaster_file}" \\
  --mito-blast "${mito_blast}" \\
  --silva-tax "${taxonomy_table}" \\
  --biof-file "${biof_blast}" \\
  --output-dir "." \\
  --prefix "${mitoPrefix}" \\
  --formats "${mitoFormats}" \\
  --min-pident ${mitoMinPident} \\
  --min-percov ${mitoMinPercov} \\
  --mitochondria-substring "${mitoMitoSubstring}" \\
  --feature-col "${mitoFeatureCol}" \\
  --taxon-col "${mitoTaxonCol}" \\
  --consensus-col "${mitoConsensusCol}" \\
  --steps "${mitoSteps}" \\
  --host-first-step "${mitoHostFirstStep}" \\
  --figsize "${mitoFigsize}" \\
  --style "${mitoStyle}" \\
  --dpi ${mitoDpi} \\
  --overwrite${noPlotsFlag}
"""
}
process FILTER_COUNTS {
    cpus pipelineThreads
    conda "${filterCountsCondaEnvPath}"
    publishDir dirMap.asv, mode: 'copy', pattern: '*', saveAs: { filename ->
        filename.endsWith('.mito.tsv') ? null : filename
    }
    publishDir filterCountsMitoDir, mode: 'copy', pattern: '*.mito.tsv'

    input:
    tuple path(count_table), path(asv_fasta)
    path(taxonomy_table)
    path(nontarget_table)

    output:
    path("${filterCountsOutputName}"), emit: filtered_counts
    path("${filterCountsOutputName}".replace('.tsv','.micro.tsv')), optional: true, emit: filtered_micro
    path("${filterCountsOutputName}".replace('.tsv','.mito.tsv')), emit: filtered_mito
    path("${filterCountsOutputName}".replace('.tsv','.decon.tsv')), optional: true, emit: filtered_decon

    when:
    filterCountsEnabled

    script:
    def aspireCacheGeneration = task.ext.aspire_cache_generation ?: 0
    def metadataArg = filterCountsMetadataPath ? """  --metadata "${filterCountsMetadataPath}" \\\n""" : ''
    def groupArg = filterCountsGroupCol ? """  --group-col "${filterCountsGroupCol}" \\\n""" : ''
    def saveInterArg = filterCountsSaveIntermediates ? "  --save-intermediates \\\n" : ''
    def mitoColsArg = (filterCountsMitoCols && !filterCountsMitoCols.isEmpty()) ?
        """  --mito-cols ${filterCountsMitoCols.collect { "\"${it}\"" }.join(' ')} \\\n""" : ''
    def excludeTaxaArg = (filterCountsExcludeTaxa && !filterCountsExcludeTaxa.isEmpty()) ?
        filterCountsExcludeTaxa.collect { item -> """  --exclude-taxon "${item}" \\\n""" }.join('') : ''
"""
set -euo pipefail
echo "filter_nontarget.py md5: ${filterCountsScriptHash}"
"\${CONDA_PREFIX}/bin/python" "${filterCountsScriptPath}" \\
  --count-table "${count_table}" \\
  --nontarget-table "${nontarget_table}" \\
  --taxonomy-table "${taxonomy_table}" \\
${metadataArg}${groupArg}  --min-group-size ${filterCountsMinGroup} \\
  --abundance-threshold ${filterCountsAbundance} \\
  --sample-id-col "${filterCountsSampleCol}" \\
  --min-consensus ${filterCountsMinConsensus} \\
  --taxon-col "${filterCountsTaxonCol}" \\
  --consensus-col "${filterCountsConsensusCol}" \\
  --biofactorial-col "${filterCountsBiofactorialCol}" \\
${mitoColsArg}${excludeTaxaArg}  --mito-output-dir "." \\
  --output "${filterCountsOutputName}" \\
${saveInterArg}
"""
}
