#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def projectRootDir = new File(projectDir.toString())
def defaultConfigPath = "${projectDir}/asv_pipeline_nextflow.yml"
def paramsMap = [:]

try {
    paramsMap = workflow.params ? new LinkedHashMap(workflow.params) : [:]
} catch( Throwable ignored ) {
    paramsMap = [:]
}

def inlineKeys = [
    'paths','resources','fastp','merge','filter','unoise',
    'table_filter','filename_patterns','environments','config_root',
    'pipeline_config'
]
def hasInlineConfig = inlineKeys.any { paramsMap.containsKey(it) }

def config
File configFile = null
File configRoot = projectRootDir

def aspirePhase = (System.getenv('ASPIRE_PHASE') ?: 'all').trim().toLowerCase()
if( !(aspirePhase in ['all', 'preprocess', 'analysis']) ) {
    exit 1, "ASPIRE_PHASE must be one of: all, preprocess, analysis (received: ${aspirePhase})"
}
log.info "ASPIRE execution phase: ${aspirePhase}"

def providedConfigPath = null
def providedConfigEnv = System.getenv('ASPIRE_PIPELINE_CONFIG') ?: System.getenv('SPARK_PIPELINE_CONFIG')
if( providedConfigEnv ) {
    providedConfigPath = providedConfigEnv
} else if( paramsMap.containsKey('pipeline_config') ) {
    providedConfigPath = paramsMap.pipeline_config
} else if( paramsMap.containsKey('config') ) {
    providedConfigPath = paramsMap.config
}

if( providedConfigPath ) {
    def configPath = file(providedConfigPath)
    configFile = configPath.toFile()
    if( !configFile.exists() ) {
        exit 1, "Config file not found: ${configFile}"
    }
    config = new groovy.yaml.YamlSlurper().parse(configFile)
    configRoot = configFile.parentFile ?: projectRootDir
    log.info "Loaded config from ${configFile}"
} else if( hasInlineConfig ) {
    config = paramsMap
    if( paramsMap.containsKey('config_root') ) {
        def rootPath = file(paramsMap.config_root)
        configRoot = rootPath.toFile()
    }
    log.info "Using inline Nextflow parameters as configuration."
} else {
    def configPath = file(defaultConfigPath)
    configFile = configPath.toFile()
    if( !configFile.exists() ) {
        exit 1, "Config file not found: ${defaultConfigPath}"
    }
    config = new groovy.yaml.YamlSlurper().parse(configFile)
    configRoot = configFile.parentFile ?: projectRootDir
    log.info "Loaded default config from ${configFile}"
}

def resolvePath = { String pathValue ->
    if( !pathValue ) {
        return null
    }
    def candidate = new File(pathValue)
    if( candidate.isAbsolute() ) {
        return candidate.canonicalPath
    }
    return new File(configRoot, pathValue).canonicalPath
}

def fileMd5(File inputFile) {
    return java.security.MessageDigest.getInstance('MD5').digest(inputFile.bytes).encodeHex().toString()
}

def r1Tokens = normalizeList(config.filename_patterns?.r1_tokens, ['R1','1'])
def r2Tokens = normalizeList(config.filename_patterns?.r2_tokens, ['R2','2'])
assert r1Tokens.size() == r2Tokens.size() : "R1 token count (${r1Tokens.size()}) must match R2 token count (${r2Tokens.size()})"

def extPatterns = compilePatterns(config.filename_patterns?.ext_patterns, ['\\.fastq\\.gz$','\\.fq\\.gz$','\\.fastq$','\\.fq$'])
def stripRegex = config.filename_patterns?.sample_strip_regex ?: '(_S[0-9]+)?(_L[0-9]{3})?(_R[12])?(_[12])?(_001)?$'

def inputDir = resolvePath(config.paths?.input_dir)
def publicOutputDir = resolvePath(config.paths?.output_dir)
assert inputDir : "paths.input_dir must be provided in the YAML config"
assert publicOutputDir : "paths.output_dir must be provided in the YAML config"
def runtimeDir = config.paths?.runtime_dir ? resolvePath(config.paths.runtime_dir) : new File(publicOutputDir, '.aspire').canonicalPath
def outputDir = new File(runtimeDir, 'publication_staging').canonicalPath
def workDirPath = config.paths?.work_dir ? resolvePath(config.paths.work_dir) : new File(runtimeDir, 'nf_work').canonicalPath
def resolvedCondaCacheDir = config.paths?.conda_cache_dir ? resolvePath(config.paths.conda_cache_dir) : new File(runtimeDir, 'conda_cache').canonicalPath
def workDirFile = new File(workDirPath)
workDirFile.mkdirs()
workflow.workDir = java.nio.file.Paths.get(workDirFile.canonicalPath)
log.info "Using Nextflow work directory: ${workflow.workDir}"
log.info "Using publication staging directory: ${outputDir}"
def condaCacheDirFile = new File(resolvedCondaCacheDir)
condaCacheDirFile.mkdirs()
System.setProperty('NXF_CONDA_CACHEDIR', condaCacheDirFile.canonicalPath)
log.info "Using custom Conda cache directory: ${condaCacheDirFile.canonicalPath}"

def allowSingleEnd = (config.resources?.single_end ?: false) as boolean
int hostThreads = Runtime.runtime.availableProcessors()
int sampleThreads = config.resources?.sample_threads ? (config.resources.sample_threads as int) :
    (config.resources?.threads ? (config.resources.threads as int) : 1)
int pipelineThreads = config.resources?.analysis_threads ? (config.resources.analysis_threads as int) : hostThreads
int maxParallelSampleTasks = config.resources?.max_parallel_sample_tasks ?
    (config.resources.max_parallel_sample_tasks as int) : Math.max(1, (pipelineThreads / sampleThreads) as int)
if( sampleThreads < 1 || pipelineThreads < 1 || maxParallelSampleTasks < 1 ) {
    exit 1, "resources.sample_threads, resources.analysis_threads, and resources.max_parallel_sample_tasks must be positive integers"
}

if( config.fastp && !(config.fastp instanceof Map) ) {
    log.warn "Ignoring non-map fastp configuration (${config.fastp.getClass()?.simpleName})"
}
def fastpConfigMap = (config.fastp instanceof Map) ? config.fastp : [:]
def fastpTrimValues = [
    front_r1: fastpConfigMap.trim_front_r1 != null ? (fastpConfigMap.trim_front_r1 as int) : 0,
    tail_r1 : fastpConfigMap.trim_tail_r1  != null ? (fastpConfigMap.trim_tail_r1  as int) : 0,
    front_r2: fastpConfigMap.trim_front_r2 != null ? (fastpConfigMap.trim_front_r2 as int) : 0,
    tail_r2 : fastpConfigMap.trim_tail_r2  != null ? (fastpConfigMap.trim_tail_r2  as int) : 0
]
if( config.merge && !(config.merge instanceof Map) ) {
    log.warn "Ignoring non-map merge configuration (${config.merge.getClass()?.simpleName})"
}
def mergeConfigMap = (config.merge instanceof Map) ? config.merge : [:]
def mergeMaxDiffs = mergeConfigMap.max_diffs != null ? (mergeConfigMap.max_diffs as int) : 20
def mergeMinOverlap = mergeConfigMap.min_overlap != null ? (mergeConfigMap.min_overlap as int) : 5
def mergeTruncQuality = mergeConfigMap.trunc_quality != null ? (mergeConfigMap.trunc_quality as int) : 5
boolean mergeAllowStagger = (mergeConfigMap.allow_stagger ?: false) as boolean

def dirMap = [
    fastp    : "${outputDir}/fastp",
    merge    : "${outputDir}/merged",
    filter   : "${outputDir}/filtered",
    concat   : "${outputDir}/concat",
    derep    : "${outputDir}/derep",
    sina     : "${outputDir}/sina",
    denoise  : "${outputDir}/denoise",
    nochi    : "${outputDir}/nochimeras",
    asv      : "${outputDir}/ASVs",
    mito     : "${outputDir}/mito",
    taxonomy : "${outputDir}/taxonomy",
    stats    : "${outputDir}/stats",
    logs     : "${outputDir}/logs",
    metadata : "${outputDir}/metadata",
    reference: "${outputDir}/reference"
]
new File(dirMap.concat).mkdirs()
new File(dirMap.sina).mkdirs()
new File(dirMap.asv).mkdirs()
new File(dirMap.mito).mkdirs()
new File(dirMap.taxonomy).mkdirs()
new File(dirMap.stats).mkdirs()
new File(dirMap.logs).mkdirs()
new File(dirMap.metadata).mkdirs()
new File(dirMap.reference).mkdirs()

def sinaReferenceFilename = 'SILVA_138.2_SSURef_NR99_03_07_24_opt.arb'
def defaultSinaReferenceUrl = 'https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_03_07_24_opt.arb.gz'
def defaultTaxonomyTaxUrl = 'https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza'
def defaultTaxonomySeqsUrl = 'https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza'
def defaultTaxonomyTaxFilename = 'silva-138_2-ssu-nr99-tax.qza'
def defaultTaxonomySeqsFilename = 'silva-138_2-ssu-nr99-seqs-DNA.qza'

def envConfigPath = config.environments?.main
def resolvedEnvPath = envConfigPath ? resolveOptionalPath(envConfigPath, configRoot) : null
def defaultEnvPath = new File("${projectDir}/processes/shared_envs/asv_pipeline.yml").canonicalPath
def condaEnvPath = resolvedEnvPath ?: defaultEnvPath
def condaEnvFile = file(condaEnvPath)
if( !condaEnvFile.exists() ) {
    exit 1, "Conda environment YAML not found: ${condaEnvPath}"
}
log.info "Using Conda/Mamba env definition: ${condaEnvPath}"

def sinaEnvConfigPath = config.environments?.sina
def resolvedSinaEnvPath = sinaEnvConfigPath ? resolveOptionalPath(sinaEnvConfigPath, configRoot) : null
def defaultSinaEnvPath = new File("${projectDir}/processes/shared_envs/sina.yml").canonicalPath
def sinaCondaEnvPath = resolvedSinaEnvPath ?: defaultSinaEnvPath
def sinaEnvFile = file(sinaCondaEnvPath)
if( !sinaEnvFile.exists() ) {
    exit 1, "SINA conda environment YAML not found: ${sinaCondaEnvPath}"
}
log.info "Using SINA Conda/Mamba env definition: ${sinaCondaEnvPath}"

def taxonomyEnvConfigPath = config.environments?.taxonomy
def resolvedTaxonomyEnvPath = taxonomyEnvConfigPath ? resolveOptionalPath(taxonomyEnvConfigPath, configRoot) : null
def defaultTaxonomyEnvPath = new File("${projectDir}/processes/shared_envs/qiime2.yml").canonicalPath
def taxonomyCondaEnvPath = resolvedTaxonomyEnvPath ?: defaultTaxonomyEnvPath
def taxonomyEnvFile = file(taxonomyCondaEnvPath)
if( !taxonomyEnvFile.exists() ) {
    exit 1, "Taxonomy conda environment YAML not found: ${taxonomyCondaEnvPath}"
}
log.info "Using taxonomy Conda/Mamba env definition: ${taxonomyCondaEnvPath}"

def mitomasterEnvConfigPath = config.environments?.mitomaster
def resolvedMitomasterEnvPath = mitomasterEnvConfigPath ? resolveOptionalPath(mitomasterEnvConfigPath, configRoot) : null
def defaultMitomasterEnvPath = new File("${projectDir}/processes/shared_envs/mitomaster.yml").canonicalPath
def mitomasterCondaEnvPath = resolvedMitomasterEnvPath ?: defaultMitomasterEnvPath
def mitomasterEnvFile = file(mitomasterCondaEnvPath)
if( !mitomasterEnvFile.exists() ) {
    exit 1, "MITOMASTER conda environment YAML not found: ${mitomasterCondaEnvPath}"
}
log.info "Using MITOMASTER Conda/Mamba env definition: ${mitomasterCondaEnvPath}"

def mitoCheckerEnvConfigPath = config.environments?.mito_checker
def resolvedMitoCheckerEnvPath = mitoCheckerEnvConfigPath ? resolveOptionalPath(mitoCheckerEnvConfigPath, configRoot) : null
def defaultMitoCheckerEnvPath = new File("${projectDir}/processes/shared_envs/mito_checker.yml").canonicalPath
def mitoCheckerCondaEnvPath = resolvedMitoCheckerEnvPath ?: defaultMitoCheckerEnvPath
def mitoCheckerEnvFile = file(mitoCheckerCondaEnvPath)
if( !mitoCheckerEnvFile.exists() ) {
    exit 1, "Mito checker conda environment YAML not found: ${mitoCheckerCondaEnvPath}"
}
log.info "Using mito checker Conda/Mamba env definition: ${mitoCheckerCondaEnvPath}"

def filterCountsEnvConfigPath = config.environments?.filter_counts
def resolvedFilterCountsEnvPath = filterCountsEnvConfigPath ? resolveOptionalPath(filterCountsEnvConfigPath, configRoot) : null
def defaultFilterCountsEnvPath = new File("${projectDir}/processes/shared_envs/filter_counts.yml").canonicalPath
def filterCountsCondaEnvPath = resolvedFilterCountsEnvPath ?: defaultFilterCountsEnvPath
def filterCountsEnvFile = file(filterCountsCondaEnvPath)
if( !filterCountsEnvFile.exists() ) {
    exit 1, "filter_counts conda environment YAML not found: ${filterCountsCondaEnvPath}"
}
log.info "Using filter_counts Conda/Mamba env definition: ${filterCountsCondaEnvPath}"

def generalStatsEnvConfigPath = config.environments?.general_stats
def resolvedGeneralStatsEnvPath = generalStatsEnvConfigPath ? resolveOptionalPath(generalStatsEnvConfigPath, configRoot) : null
def defaultGeneralStatsEnvPath = new File("${projectDir}/processes/shared_envs/general_stats.yml").canonicalPath
def generalStatsCondaEnvPath = resolvedGeneralStatsEnvPath ?: defaultGeneralStatsEnvPath
def generalStatsEnvFile = file(generalStatsCondaEnvPath)
if( !generalStatsEnvFile.exists() ) {
    exit 1, "general_stats conda environment YAML not found: ${generalStatsCondaEnvPath}"
}
log.info "Using general_stats Conda/Mamba env definition: ${generalStatsCondaEnvPath}"

def sankeyEnvConfigPath = config.environments?.sankey
def resolvedSankeyEnvPath = sankeyEnvConfigPath ? resolveOptionalPath(sankeyEnvConfigPath, configRoot) : null
def defaultSankeyEnvPath = new File("${projectDir}/processes/shared_envs/sankey.yml").canonicalPath
def sankeyCondaEnvPath = resolvedSankeyEnvPath ?: defaultSankeyEnvPath
def sankeyEnvFile = file(sankeyCondaEnvPath)
if( !sankeyEnvFile.exists() ) {
    exit 1, "Sankey conda environment YAML not found: ${sankeyCondaEnvPath}"
}
log.info "Using Sankey Conda/Mamba env definition: ${sankeyCondaEnvPath}"

def plotMetadataEnvConfigPath = config.environments?.plot_metadata
def resolvedPlotMetadataEnvPath = plotMetadataEnvConfigPath ? resolveOptionalPath(plotMetadataEnvConfigPath, configRoot) : null
def defaultPlotMetadataEnvPath = new File("${projectDir}/processes/shared_envs/plot_metadata.yml").canonicalPath
def plotMetadataCondaEnvPath = resolvedPlotMetadataEnvPath ?: defaultPlotMetadataEnvPath
def plotMetadataEnvFile = file(plotMetadataCondaEnvPath)
if( !plotMetadataEnvFile.exists() ) {
    exit 1, "Plot metadata conda environment YAML not found: ${plotMetadataCondaEnvPath}"
}
log.info "Using plot metadata Conda/Mamba env definition: ${plotMetadataCondaEnvPath}"

def asvTimeDepthCurtainEnvConfigPath = config.environments?.asv_time_depth_curtain
def resolvedAsvTimeDepthCurtainEnvPath = asvTimeDepthCurtainEnvConfigPath ? resolveOptionalPath(asvTimeDepthCurtainEnvConfigPath, configRoot) : null
def defaultAsvTimeDepthCurtainEnvPath = new File("${projectDir}/processes/asv_time_depth_curtain/env.yml").canonicalPath
def asvTimeDepthCurtainCondaEnvPath = resolvedAsvTimeDepthCurtainEnvPath ?: defaultAsvTimeDepthCurtainEnvPath
if( !file(asvTimeDepthCurtainCondaEnvPath).exists() ) {
    exit 1, "ASV time-depth curtain conda environment YAML not found: ${asvTimeDepthCurtainCondaEnvPath}"
}

def batchCorrectionEnvConfigPath = config.environments?.batch_correction
def resolvedBatchCorrectionEnvPath = batchCorrectionEnvConfigPath ? resolveOptionalPath(batchCorrectionEnvConfigPath, configRoot) : null
def defaultBatchCorrectionEnvPath = new File("${projectDir}/processes/shared_envs/asv_batch_correction.yml").canonicalPath
def batchCorrectionCondaEnvPath = resolvedBatchCorrectionEnvPath ?: defaultBatchCorrectionEnvPath
def batchCorrectionEnvFile = file(batchCorrectionCondaEnvPath)
if( !batchCorrectionEnvFile.exists() ) {
    exit 1, "Batch correction conda environment YAML not found: ${batchCorrectionCondaEnvPath}"
}
log.info "Using batch correction Conda/Mamba env definition: ${batchCorrectionCondaEnvPath}"

def outlierEnvConfigPath = config.environments?.outlier_checker
def resolvedOutlierEnvPath = outlierEnvConfigPath ? resolveOptionalPath(outlierEnvConfigPath, configRoot) : null
def defaultOutlierEnvPath = new File("${projectDir}/processes/shared_envs/outlier_checker.yml").canonicalPath
def outlierCondaEnvPath = resolvedOutlierEnvPath ?: defaultOutlierEnvPath
def outlierEnvFile = file(outlierCondaEnvPath)
if( !outlierEnvFile.exists() ) {
    exit 1, "Outlier checker conda environment YAML not found: ${outlierCondaEnvPath}"
}
log.info "Using outlier checker Conda/Mamba env definition: ${outlierCondaEnvPath}"

def collectorsEnvConfigPath = config.environments?.collectors_curve
def resolvedCollectorsEnvPath = collectorsEnvConfigPath ? resolveOptionalPath(collectorsEnvConfigPath, configRoot) : null
def defaultCollectorsEnvPath = new File("${projectDir}/processes/shared_envs/collectors_curve.yml").canonicalPath
def collectorsCondaEnvPath = resolvedCollectorsEnvPath ?: defaultCollectorsEnvPath
def collectorsEnvFile = file(collectorsCondaEnvPath)
if( !collectorsEnvFile.exists() ) {
    exit 1, "Collectors curve conda environment YAML not found: ${collectorsCondaEnvPath}"
}
log.info "Using collectors curve Conda/Mamba env definition: ${collectorsCondaEnvPath}"

def diversityEnvConfigPath = config.environments?.diversity
def resolvedDiversityEnvPath = diversityEnvConfigPath ? resolveOptionalPath(diversityEnvConfigPath, configRoot) : null
def diversityCondaEnvPath = resolvedDiversityEnvPath ?: new File("${projectDir}/processes/diversity_analysis/env.yml").canonicalPath
def diversityEnvFile = file(diversityCondaEnvPath)
if( !diversityEnvFile.exists() ) {
    exit 1, "Diversity conda environment YAML not found: ${diversityCondaEnvPath}"
}
log.info "Using diversity Conda/Mamba env definition: ${diversityCondaEnvPath}"

def indicspeciesEnvConfigPath = config.environments?.indicspecies
def resolvedIndicspeciesEnvPath = indicspeciesEnvConfigPath ? resolveOptionalPath(indicspeciesEnvConfigPath, configRoot) : null
def indicspeciesCondaEnvPath = resolvedIndicspeciesEnvPath ?: new File("${projectDir}/processes/indicspecies/env.yml").canonicalPath
def indicspeciesEnvFile = file(indicspeciesCondaEnvPath)
if( !indicspeciesEnvFile.exists() ) {
    exit 1, "Indicspecies conda environment YAML not found: ${indicspeciesCondaEnvPath}"
}
log.info "Using indicspecies Conda/Mamba env definition: ${indicspeciesCondaEnvPath}"

def clustermapsEnvConfigPath = config.environments?.clustermaps
def resolvedClustermapsEnvPath = clustermapsEnvConfigPath ? resolveOptionalPath(clustermapsEnvConfigPath, configRoot) : null
def clustermapsCondaEnvPath = resolvedClustermapsEnvPath ?: new File("${projectDir}/processes/clustermaps/env.yml").canonicalPath
def clustermapsEnvFile = file(clustermapsCondaEnvPath)
if( !clustermapsEnvFile.exists() ) {
    exit 1, "Clustermaps conda environment YAML not found: ${clustermapsCondaEnvPath}"
}
log.info "Using clustermaps Conda/Mamba env definition: ${clustermapsCondaEnvPath}"

def spieceasiEnvConfigPath = config.environments?.spieceasi
def resolvedSpieceasiEnvPath = spieceasiEnvConfigPath ? resolveOptionalPath(spieceasiEnvConfigPath, configRoot) : null
def spieceasiCondaEnvPath = resolvedSpieceasiEnvPath ?: new File("${projectDir}/processes/spieceasi/env.yml").canonicalPath
def spieceasiEnvFile = file(spieceasiCondaEnvPath)
if( !spieceasiEnvFile.exists() ) {
    exit 1, "SPIEC-EASI conda environment YAML not found: ${spieceasiCondaEnvPath}"
}
log.info "Using SPIEC-EASI Conda/Mamba env definition: ${spieceasiCondaEnvPath}"

def networkEnvConfigPath = config.environments?.network
def resolvedNetworkEnvPath = networkEnvConfigPath ? resolveOptionalPath(networkEnvConfigPath, configRoot) : null
def networkCondaEnvPath = resolvedNetworkEnvPath ?: new File("${projectDir}/processes/graph_network/env.yml").canonicalPath
def networkEnvFile = file(networkCondaEnvPath)
if( !networkEnvFile.exists() ) {
    exit 1, "Network conda environment YAML not found: ${networkCondaEnvPath}"
}
log.info "Using network Conda/Mamba env definition: ${networkCondaEnvPath}"

def networkModulesEnvConfigPath = config.environments?.network_modules
def resolvedNetworkModulesEnvPath = networkModulesEnvConfigPath ? resolveOptionalPath(networkModulesEnvConfigPath, configRoot) : null
def networkModulesCondaEnvPath = resolvedNetworkModulesEnvPath ?: spieceasiCondaEnvPath
def networkModulesEnvFile = file(networkModulesCondaEnvPath)
if( !networkModulesEnvFile.exists() ) {
    exit 1, "Network modules conda environment YAML not found: ${networkModulesCondaEnvPath}"
}
log.info "Using network modules Conda/Mamba env definition: ${networkModulesCondaEnvPath}"

def genomeCooccurrenceEnvConfigPath = config.environments?.genome_cooccurrence
def resolvedGenomeCooccurrenceEnvPath = genomeCooccurrenceEnvConfigPath ? resolveOptionalPath(genomeCooccurrenceEnvConfigPath, configRoot) : null
def genomeCooccurrenceCondaEnvPath = resolvedGenomeCooccurrenceEnvPath ?: new File("${projectDir}/processes/genome_cooccurrence/env.yml").canonicalPath
def genomeCooccurrenceEnvFile = file(genomeCooccurrenceCondaEnvPath)
if( !genomeCooccurrenceEnvFile.exists() ) {
    exit 1, "Genome co-occurrence conda environment YAML not found: ${genomeCooccurrenceCondaEnvPath}"
}
log.info "Using genome co-occurrence Conda/Mamba env definition: ${genomeCooccurrenceCondaEnvPath}"

def masterSummaryEnvConfigPath = config.environments?.master_summary
def resolvedMasterSummaryEnvPath = masterSummaryEnvConfigPath ? resolveOptionalPath(masterSummaryEnvConfigPath, configRoot) : null
def masterSummaryCondaEnvPath = resolvedMasterSummaryEnvPath ?: new File("${projectDir}/processes/master_summary/env.yml").canonicalPath
def masterSummaryEnvFile = file(masterSummaryCondaEnvPath)
if( !masterSummaryEnvFile.exists() ) {
    exit 1, "Master summary conda environment YAML not found: ${masterSummaryCondaEnvPath}"
}
log.info "Using master summary Conda/Mamba env definition: ${masterSummaryCondaEnvPath}"

def asvMagLinkEnvConfigPath = config.environments?.asv_mag_link
def resolvedAsvMagLinkEnvPath = asvMagLinkEnvConfigPath ? resolveOptionalPath(asvMagLinkEnvConfigPath, configRoot) : null
def asvMagLinkCondaEnvPath = resolvedAsvMagLinkEnvPath ?: new File("${projectDir}/processes/asv_mag_link/env.yml").canonicalPath
def asvMagLinkEnvFile = file(asvMagLinkCondaEnvPath)
if( !asvMagLinkEnvFile.exists() ) {
    exit 1, "ASV-MAG linking conda environment YAML not found: ${asvMagLinkCondaEnvPath}"
}
log.info "Using ASV-MAG linking Conda/Mamba env definition: ${asvMagLinkCondaEnvPath}"

def asvMagNetworkEnvConfigPath = config.environments?.asv_mag_network
def resolvedAsvMagNetworkEnvPath = asvMagNetworkEnvConfigPath ? resolveOptionalPath(asvMagNetworkEnvConfigPath, configRoot) : null
def asvMagNetworkCondaEnvPath = resolvedAsvMagNetworkEnvPath ?: new File("${projectDir}/processes/asv_mag_network/env.yml").canonicalPath
def asvMagNetworkEnvFile = file(asvMagNetworkCondaEnvPath)
if( !asvMagNetworkEnvFile.exists() ) {
    exit 1, "ASV-MAG network conda environment YAML not found: ${asvMagNetworkCondaEnvPath}"
}
log.info "Using ASV-MAG network Conda/Mamba env definition: ${asvMagNetworkCondaEnvPath}"

def asvMagCurtainsEnvConfigPath = config.environments?.asv_mag_curtains
def resolvedAsvMagCurtainsEnvPath = asvMagCurtainsEnvConfigPath ? resolveOptionalPath(asvMagCurtainsEnvConfigPath, configRoot) : null
def asvMagCurtainsCondaEnvPath = resolvedAsvMagCurtainsEnvPath ?: new File("${projectDir}/processes/asv_mag_curtains/env.yml").canonicalPath
if( !file(asvMagCurtainsCondaEnvPath).exists() ) {
    exit 1, "ASV-MAG curtain Conda/Mamba environment YAML not found: ${asvMagCurtainsCondaEnvPath}"
}
log.info "Using ASV-MAG curtain Conda/Mamba env definition: ${asvMagCurtainsCondaEnvPath}"

def groupGuildFunctionEnvConfigPath = config.environments?.group_guild_function
def resolvedGroupGuildFunctionEnvPath = groupGuildFunctionEnvConfigPath ? resolveOptionalPath(groupGuildFunctionEnvConfigPath, configRoot) : null
def groupGuildFunctionCondaEnvPath = resolvedGroupGuildFunctionEnvPath ?: new File("${projectDir}/processes/group_guild_function/env.yml").canonicalPath
if( !file(groupGuildFunctionCondaEnvPath).exists() ) {
    exit 1, "Group-guild-function environment YAML not found: ${groupGuildFunctionCondaEnvPath}"
}

def powerAnalysisEnvConfigPath = config.environments?.group_power_analysis ?: config.environments?.power_analysis
def resolvedPowerAnalysisEnvPath = powerAnalysisEnvConfigPath ? resolveOptionalPath(powerAnalysisEnvConfigPath, configRoot) : null
def powerAnalysisCondaEnvPath = resolvedPowerAnalysisEnvPath ?: new File("${projectDir}/processes/power_analysis_pipeline/env.yml").canonicalPath
def powerAnalysisEnvFile = file(powerAnalysisCondaEnvPath)
if( !powerAnalysisEnvFile.exists() ) {
    exit 1, "Power analysis conda environment YAML not found: ${powerAnalysisCondaEnvPath}"
}
log.info "Using power analysis Conda/Mamba env definition: ${powerAnalysisCondaEnvPath}"

def taxonomyPatientAwareEnvConfigPath = config.environments?.taxonomy_group_association ?: config.environments?.taxonomy_patient_aware
def resolvedTaxonomyPatientAwareEnvPath = taxonomyPatientAwareEnvConfigPath ? resolveOptionalPath(taxonomyPatientAwareEnvConfigPath, configRoot) : null
def taxonomyPatientAwareCondaEnvPath = resolvedTaxonomyPatientAwareEnvPath ?: new File("${projectDir}/processes/taxonomy_patient_aware/env.yml").canonicalPath
def taxonomyPatientAwareEnvFile = file(taxonomyPatientAwareCondaEnvPath)
if( !taxonomyPatientAwareEnvFile.exists() ) {
    exit 1, "Taxonomy patient-aware conda environment YAML not found: ${taxonomyPatientAwareCondaEnvPath}"
}
log.info "Using taxonomy patient-aware Conda/Mamba env definition: ${taxonomyPatientAwareCondaEnvPath}"

def lungStatusAnalysisEnvConfigPath = config.environments?.paired_group_contrast ?: config.environments?.lung_status_analysis
def resolvedLungStatusAnalysisEnvPath = lungStatusAnalysisEnvConfigPath ? resolveOptionalPath(lungStatusAnalysisEnvConfigPath, configRoot) : null
def lungStatusAnalysisCondaEnvPath = resolvedLungStatusAnalysisEnvPath ?: new File("${projectDir}/processes/lung_status_analysis/env.yml").canonicalPath
def lungStatusAnalysisEnvFile = file(lungStatusAnalysisCondaEnvPath)
if( !lungStatusAnalysisEnvFile.exists() ) {
    exit 1, "Lung status analysis conda environment YAML not found: ${lungStatusAnalysisCondaEnvPath}"
}
log.info "Using lung status analysis Conda/Mamba env definition: ${lungStatusAnalysisCondaEnvPath}"

def plotUpsetEnvConfigPath = config.environments?.plot_upset
def resolvedPlotUpsetEnvPath = plotUpsetEnvConfigPath ? resolveOptionalPath(plotUpsetEnvConfigPath, configRoot) : null
def plotUpsetCondaEnvPath = resolvedPlotUpsetEnvPath ?: new File("${projectDir}/processes/plot_upset/env.yml").canonicalPath
def plotUpsetEnvFile = file(plotUpsetCondaEnvPath)
if( !plotUpsetEnvFile.exists() ) {
    exit 1, "Plot Upset conda environment YAML not found: ${plotUpsetCondaEnvPath}"
}
log.info "Using Plot Upset Conda/Mamba env definition: ${plotUpsetCondaEnvPath}"

def bubbleplotterEnvConfigPath = config.environments?.bubbleplotter
def resolvedBubbleplotterEnvPath = bubbleplotterEnvConfigPath ? resolveOptionalPath(bubbleplotterEnvConfigPath, configRoot) : null
def bubbleplotterCondaEnvPath = resolvedBubbleplotterEnvPath ?: new File("${projectDir}/processes/bubbleplotter/env.yml").canonicalPath
def bubbleplotterEnvFile = file(bubbleplotterCondaEnvPath)
if( !bubbleplotterEnvFile.exists() ) {
    exit 1, "Bubbleplotter conda environment YAML not found: ${bubbleplotterCondaEnvPath}"
}
log.info "Using bubbleplotter Conda/Mamba env definition: ${bubbleplotterCondaEnvPath}"

def umapClusteringEnvConfigPath = config.environments?.umap_clustering
def resolvedUmapClusteringEnvPath = umapClusteringEnvConfigPath ? resolveOptionalPath(umapClusteringEnvConfigPath, configRoot) : null
def umapClusteringCondaEnvPath = resolvedUmapClusteringEnvPath ?: new File("${projectDir}/processes/umap_clustering/env.yml").canonicalPath
def umapClusteringEnvFile = file(umapClusteringCondaEnvPath)
if( !umapClusteringEnvFile.exists() ) {
    exit 1, "UMAP clustering conda environment YAML not found: ${umapClusteringCondaEnvPath}"
}
log.info "Using UMAP clustering Conda/Mamba env definition: ${umapClusteringCondaEnvPath}"

def vocCorrelationEnvConfigPath = config.environments?.voc_correlation
def resolvedVocCorrelationEnvPath = vocCorrelationEnvConfigPath ? resolveOptionalPath(vocCorrelationEnvConfigPath, configRoot) : null
def defaultVocCorrelationEnvPath = new File("${projectDir}/processes/voc_correlation/env.yml").canonicalPath
def vocCorrelationCondaEnvPath = resolvedVocCorrelationEnvPath ?: defaultVocCorrelationEnvPath
def vocCorrelationEnvFile = file(vocCorrelationCondaEnvPath)
if( !vocCorrelationEnvFile.exists() ) {
    exit 1, "VOC correlation conda environment YAML not found: ${vocCorrelationCondaEnvPath}"
}
log.info "Using VOC correlation Conda/Mamba env definition: ${vocCorrelationCondaEnvPath}"

def measurementAssociationEnvConfigPath = config.environments?.measurement_association
def resolvedMeasurementAssociationEnvPath = measurementAssociationEnvConfigPath ? resolveOptionalPath(measurementAssociationEnvConfigPath, configRoot) : null
def defaultMeasurementAssociationEnvPath = new File("${projectDir}/processes/measurement_association/env.yml").canonicalPath
def measurementAssociationCondaEnvPath = resolvedMeasurementAssociationEnvPath ?: defaultMeasurementAssociationEnvPath
def measurementAssociationEnvFile = file(measurementAssociationCondaEnvPath)
if( !measurementAssociationEnvFile.exists() ) {
    exit 1, "Measurement association conda environment YAML not found: ${measurementAssociationCondaEnvPath}"
}
log.info "Using measurement association Conda/Mamba env definition: ${measurementAssociationCondaEnvPath}"

def titanEnvConfigPath = config.environments?.titan
def resolvedTitanEnvPath = titanEnvConfigPath ? resolveOptionalPath(titanEnvConfigPath, configRoot) : null
def defaultTitanEnvPath = new File("${projectDir}/processes/titan/env.yml").canonicalPath
def titanCondaEnvPath = resolvedTitanEnvPath ?: defaultTitanEnvPath
if( !file(titanCondaEnvPath).exists() ) {
    exit 1, "TITAN Conda/Mamba environment YAML not found: ${titanCondaEnvPath}"
}
log.info "Using TITAN Conda/Mamba env definition: ${titanCondaEnvPath}"

def microbialCompartmentEnvConfigPath = config.environments?.microbial_compartments
def resolvedMicrobialCompartmentEnvPath = microbialCompartmentEnvConfigPath ? resolveOptionalPath(microbialCompartmentEnvConfigPath, configRoot) : null
def defaultMicrobialCompartmentEnvPath = new File("${projectDir}/processes/microbial_compartments/env.yml").canonicalPath
def microbialCompartmentCondaEnvPath = resolvedMicrobialCompartmentEnvPath ?: defaultMicrobialCompartmentEnvPath
if( !file(microbialCompartmentCondaEnvPath).exists() ) {
    exit 1, "Microbial-compartment Conda/Mamba environment YAML not found: ${microbialCompartmentCondaEnvPath}"
}
log.info "Using microbial-compartment Conda/Mamba env definition: ${microbialCompartmentCondaEnvPath}"

def microbialStateInterpretationEnvConfigPath = config.environments?.microbial_state_interpretation
def resolvedMicrobialStateInterpretationEnvPath = microbialStateInterpretationEnvConfigPath ? resolveOptionalPath(microbialStateInterpretationEnvConfigPath, configRoot) : null
def defaultMicrobialStateInterpretationEnvPath = new File("${projectDir}/processes/microbial_state_interpretation/env.yml").canonicalPath
def microbialStateInterpretationCondaEnvPath = resolvedMicrobialStateInterpretationEnvPath ?: defaultMicrobialStateInterpretationEnvPath
if( !file(microbialStateInterpretationCondaEnvPath).exists() ) {
    exit 1, "Microbial-state interpretation Conda/Mamba environment YAML not found: ${microbialStateInterpretationCondaEnvPath}"
}
log.info "Using microbial-state interpretation Conda/Mamba env definition: ${microbialStateInterpretationCondaEnvPath}"

def ecologicalContextAtlasEnvConfigPath = config.environments?.ecological_context_atlas
def resolvedEcologicalContextAtlasEnvPath = ecologicalContextAtlasEnvConfigPath ? resolveOptionalPath(ecologicalContextAtlasEnvConfigPath, configRoot) : null
def defaultEcologicalContextAtlasEnvPath = new File("${projectDir}/processes/ecological_context_atlas/env.yml").canonicalPath
def ecologicalContextAtlasCondaEnvPath = resolvedEcologicalContextAtlasEnvPath ?: defaultEcologicalContextAtlasEnvPath
if( !file(ecologicalContextAtlasCondaEnvPath).exists() ) {
    exit 1, "Ecological-context atlas Conda/Mamba environment YAML not found: ${ecologicalContextAtlasCondaEnvPath}"
}
log.info "Using ecological-context atlas Conda/Mamba env definition: ${ecologicalContextAtlasCondaEnvPath}"

def communityTurnoverEnvConfigPath = config.environments?.community_turnover
def resolvedCommunityTurnoverEnvPath = communityTurnoverEnvConfigPath ? resolveOptionalPath(communityTurnoverEnvConfigPath, configRoot) : null
def defaultCommunityTurnoverEnvPath = new File("${projectDir}/processes/community_turnover/env.yml").canonicalPath
def communityTurnoverCondaEnvPath = resolvedCommunityTurnoverEnvPath ?: defaultCommunityTurnoverEnvPath
if( !file(communityTurnoverCondaEnvPath).exists() ) {
    exit 1, "Community-turnover Conda/Mamba environment YAML not found: ${communityTurnoverCondaEnvPath}"
}
log.info "Using community-turnover Conda/Mamba env definition: ${communityTurnoverCondaEnvPath}"

def groupingDiagnosticsEnvConfigPath = config.environments?.grouping_diagnostics
def resolvedGroupingDiagnosticsEnvPath = groupingDiagnosticsEnvConfigPath ? resolveOptionalPath(groupingDiagnosticsEnvConfigPath, configRoot) : null
def defaultGroupingDiagnosticsEnvPath = new File("${projectDir}/processes/grouping_diagnostics/env.yml").canonicalPath
def groupingDiagnosticsCondaEnvPath = resolvedGroupingDiagnosticsEnvPath ?: defaultGroupingDiagnosticsEnvPath
def groupingDiagnosticsEnvFile = file(groupingDiagnosticsCondaEnvPath)
if( !groupingDiagnosticsEnvFile.exists() ) {
    exit 1, "Grouping diagnostics conda environment YAML not found: ${groupingDiagnosticsCondaEnvPath}"
}
log.info "Using grouping diagnostics Conda/Mamba env definition: ${groupingDiagnosticsCondaEnvPath}"

def communityPredictorEnvConfigPath = config.environments?.community_predictor_comparison
def resolvedCommunityPredictorEnvPath = communityPredictorEnvConfigPath ? resolveOptionalPath(communityPredictorEnvConfigPath, configRoot) : null
def defaultCommunityPredictorEnvPath = new File("${projectDir}/processes/community_predictor_comparison/env.yml").canonicalPath
def communityPredictorCondaEnvPath = resolvedCommunityPredictorEnvPath ?: defaultCommunityPredictorEnvPath
def communityPredictorEnvFile = file(communityPredictorCondaEnvPath)
if( !communityPredictorEnvFile.exists() ) {
    exit 1, "Community predictor comparison conda environment YAML not found: ${communityPredictorCondaEnvPath}"
}
log.info "Using community predictor comparison Conda/Mamba env definition: ${communityPredictorCondaEnvPath}"

def groupLabelAugmentationEnvConfigPath = config.environments?.group_label_augmentation
def resolvedGroupLabelAugmentationEnvPath = groupLabelAugmentationEnvConfigPath ? resolveOptionalPath(groupLabelAugmentationEnvConfigPath, configRoot) : null
def defaultGroupLabelAugmentationEnvPath = new File("${projectDir}/processes/group_label_augmentation/env.yml").canonicalPath
def groupLabelAugmentationCondaEnvPath = resolvedGroupLabelAugmentationEnvPath ?: defaultGroupLabelAugmentationEnvPath
def groupLabelAugmentationEnvFile = file(groupLabelAugmentationCondaEnvPath)
if( !groupLabelAugmentationEnvFile.exists() ) {
    exit 1, "Group-label augmentation conda environment YAML not found: ${groupLabelAugmentationCondaEnvPath}"
}
log.info "Using group-label augmentation Conda/Mamba env definition: ${groupLabelAugmentationCondaEnvPath}"

def configuredManifestPath = config.paths?.manifest ? resolveOptionalPath(config.paths.manifest, configRoot) : null
def manifestPath = configuredManifestPath
def rawReadsChannel = Channel.empty()
def sampleRecords = []
if( aspirePhase != 'analysis' ) {
    if( configuredManifestPath ) {
        sampleRecords = loadManifestSamples(configuredManifestPath)
    } else {
        sampleRecords = collectSampleRecords(inputDir, r1Tokens, r2Tokens, extPatterns, stripRegex, allowSingleEnd)
    }
    if( !sampleRecords ) {
        exit 1, configuredManifestPath ? "No usable entries detected in manifest ${configuredManifestPath}" : "No usable FASTQ files detected in ${inputDir}"
    }
    manifestPath = writeNormalizedManifest(
        sampleRecords,
        new File(dirMap.metadata, 'run_manifest.tsv'),
        configuredManifestPath
    )
    log.info "Discovered ${sampleRecords.size()} input items from ${configuredManifestPath ? "manifest ${configuredManifestPath}" : inputDir}. Normalized manifest: ${manifestPath}. Sample threads=${sampleThreads}, max parallel sample tasks=${maxParallelSampleTasks}, analysis threads=${pipelineThreads}"
    rawReadsChannel = Channel
        .from(sampleRecords)
        .map { rec ->
            def meta = [ sample_id: rec.sample_id, paired: rec.paired ]
            def r1File = file(rec.r1)
            def r2File = (rec.paired && rec.r2) ? file(rec.r2) : r1File
            tuple(meta, r1File, r2File)
        }
}

def tabFilterScript = resolveOptionalPath(config.table_filter?.script, configRoot) ?: "${projectDir}/processes/filter_table/filter_ASV_table.py"
def tableScriptFile = file(tabFilterScript)
if( !tableScriptFile.exists() ) {
    exit 1, "Table filter script not found: ${tabFilterScript}"
}

def preprocessDatasetBuildScriptFile = new File("${projectDir}/processes/preprocessing_dataset/build_preprocessing_dataset.py")
if( !preprocessDatasetBuildScriptFile.exists() ) {
    exit 1, "Canonical preprocessing dataset builder not found: ${preprocessDatasetBuildScriptFile}"
}
def preprocessDatasetBuildScriptPath = preprocessDatasetBuildScriptFile.canonicalPath
def preprocessDatasetValidateScriptFile = new File("${projectDir}/processes/preprocessing_dataset/validate_preprocessing_dataset.py")
if( !preprocessDatasetValidateScriptFile.exists() ) {
    exit 1, "Canonical preprocessing dataset validator not found: ${preprocessDatasetValidateScriptFile}"
}
def preprocessDatasetValidateScriptPath = preprocessDatasetValidateScriptFile.canonicalPath
def preprocessingParametersJson = groovy.json.JsonOutput.toJson([
    filename_patterns: config.filename_patterns ?: [:],
    resources: [single_end: allowSingleEnd, sample_threads: sampleThreads],
    fastp: config.fastp ?: [:],
    merge: config.merge ?: [:],
    filter: config.filter ?: [:],
    concat: config.concat ?: [:],
    unoise: config.unoise ?: [:],
    table_filter: config.table_filter ?: [:],
    taxonomy: config.taxonomy ?: [:],
    mito: config.mito ?: [:],
    filter_counts: config.filter_counts ?: [:]
])

def concatConfig = config.concat ?: [:]
def concatRelabelEnabled = concatConfig.containsKey('relabel') ? (concatConfig.relabel as boolean) : true
def concatLabelSep = concatConfig.label_sep ?: ':'
def filterConfig = config.filter ?: [:]
def filterMaxEe = filterConfig.max_ee != null ? (filterConfig.max_ee as double) : 1.0d
def filterMinLen = filterConfig.min_len != null ? (filterConfig.min_len as int) : 245
def filterMaxLen = filterConfig.max_len != null ? (filterConfig.max_len as int) : 1500

def parseSinaScriptFile = new File("${projectDir}/processes/sina_trim/parse_sina_log.py")
if( !parseSinaScriptFile.exists() ) {
    exit 1, "parse_sina_log.py not found in project directory"
}
def parseSinaScriptPath = parseSinaScriptFile.canonicalPath
def trimSinaScriptFile = new File("${projectDir}/processes/sina_trim/trim_v_sina.py")
if( !trimSinaScriptFile.exists() ) {
    exit 1, "trim_v_sina.py not found in project directory"
}
def trimSinaScriptPath = trimSinaScriptFile.canonicalPath
def mitomasterScriptFile = new File("${projectDir}/processes/mitomaster/mitomaster.py")
if( !mitomasterScriptFile.exists() ) {
    exit 1, "mitomaster.py not found in project directory"
}
def mitomasterScriptPath = mitomasterScriptFile.canonicalPath
def mitoCheckerScriptFile = new File("${projectDir}/processes/mito_decontam/mito_checker.py")
if( !mitoCheckerScriptFile.exists() ) {
    exit 1, "mito_checker.py not found in project directory"
}
def mitoCheckerScriptPath = mitoCheckerScriptFile.canonicalPath
def sankeyScriptFile = new File("${projectDir}/processes/sankey/sankey_builder.py")
if( !sankeyScriptFile.exists() ) {
    exit 1, "sankey_builder.py not found in project directory"
}
def sankeyScriptPath = sankeyScriptFile.canonicalPath
def sankeyScriptHash = fileMd5(sankeyScriptFile)
def plotMetadataScriptFile = new File("${projectDir}/processes/plot_metadata/plot_metadata.py")
if( !plotMetadataScriptFile.exists() ) {
    exit 1, "plot_metadata.py not found in project directory"
}
def plotMetadataScriptPath = plotMetadataScriptFile.canonicalPath
def plotMetadataScriptHash = fileMd5(plotMetadataScriptFile)
def asvTimeDepthCurtainScriptFile = new File("${projectDir}/processes/asv_time_depth_curtain/asv_time_depth_curtain.py")
if( !asvTimeDepthCurtainScriptFile.exists() ) {
    exit 1, "asv_time_depth_curtain.py not found in project directory"
}
def asvTimeDepthCurtainScriptPath = asvTimeDepthCurtainScriptFile.canonicalPath
def asvTimeDepthCurtainScriptHash = fileMd5(asvTimeDepthCurtainScriptFile)
def batchCorrectionScriptFile = new File("${projectDir}/processes/asv_batch_correction/asv_batch_correction.py")
if( !batchCorrectionScriptFile.exists() ) {
    exit 1, "asv_batch_correction.py not found in project directory"
}
def batchCorrectionScriptPath = batchCorrectionScriptFile.canonicalPath
def outlierCheckerScriptFile = new File("${projectDir}/processes/outlier_checker/outlier_checker.py")
if( !outlierCheckerScriptFile.exists() ) {
    exit 1, "outlier_checker.py not found in project directory"
}
def outlierCheckerScriptPath = outlierCheckerScriptFile.canonicalPath
def collectorsCurveScriptFile = new File("${projectDir}/processes/collectors_curve/collectors_curve.py")
if( !collectorsCurveScriptFile.exists() ) {
    exit 1, "collectors_curve.py not found in project directory"
}
def collectorsCurveScriptPath = collectorsCurveScriptFile.canonicalPath
def filterCountsScriptFile = new File("${projectDir}/processes/filter_counts/filter_nontarget.py")
if( !filterCountsScriptFile.exists() ) {
    exit 1, "filter_nontarget.py not found in project directory"
}
def filterCountsScriptPath = filterCountsScriptFile.canonicalPath
def filterCountsScriptHash = fileMd5(filterCountsScriptFile)
def calcDivScriptFile = new File("${projectDir}/processes/diversity_analysis/calc_div.py")
if( !calcDivScriptFile.exists() ) {
    exit 1, "calc_div.py not found in project directory"
}
def calcDivScriptPath = calcDivScriptFile.canonicalPath
def plotDiversityScriptFile = new File("${projectDir}/processes/diversity_analysis/plot_diversity.py")
if( !plotDiversityScriptFile.exists() ) {
    exit 1, "plot_diversity.py not found in project directory"
}
def plotDiversityScriptPath = plotDiversityScriptFile.canonicalPath
def plotUpsetScriptFile = new File("${projectDir}/processes/plot_upset/plot_upset.py")
if( !plotUpsetScriptFile.exists() ) {
    exit 1, "plot_upset.py not found in project directory"
}
def plotUpsetScriptPath = plotUpsetScriptFile.canonicalPath
def bubbleplotterScriptFile = new File("${projectDir}/processes/bubbleplotter/bubbleplotter.py")
if( !bubbleplotterScriptFile.exists() ) {
    exit 1, "bubbleplotter.py not found in project directory"
}
def bubbleplotterScriptPath = bubbleplotterScriptFile.canonicalPath
def umapClusteringScriptFile = new File("${projectDir}/processes/umap_clustering/umap_clustering.py")
if( !umapClusteringScriptFile.exists() ) {
    exit 1, "umap_clustering.py not found in project directory"
}
def umapClusteringScriptPath = umapClusteringScriptFile.canonicalPath
def indicspeciesScriptFile = new File("${projectDir}/processes/indicspecies/run_indicspecies.R")
if( !indicspeciesScriptFile.exists() ) {
    exit 1, "run_indicspecies.R not found in project directory"
}
def indicspeciesScriptPath = indicspeciesScriptFile.canonicalPath
def plotIndicspeciesScriptFile = new File("${projectDir}/processes/indicspecies_plots/plot_indicspecies.py")
if( !plotIndicspeciesScriptFile.exists() ) {
    exit 1, "plot_indicspecies.py not found in project directory"
}
def plotIndicspeciesScriptPath = plotIndicspeciesScriptFile.canonicalPath
def plotIndicspeciesAlignedScriptFile = new File("${projectDir}/processes/indicspecies_aligned_plots/plot_indicspecies_aligned.py")
if( !plotIndicspeciesAlignedScriptFile.exists() ) {
    exit 1, "plot_indicspecies_aligned.py not found in project directory"
}
def plotIndicspeciesAlignedScriptPath = plotIndicspeciesAlignedScriptFile.canonicalPath
def clustermapsScriptFile = new File("${projectDir}/processes/clustermaps/plot_clustermaps.py")
if( !clustermapsScriptFile.exists() ) {
    exit 1, "plot_clustermaps.py not found in project directory"
}
def clustermapsScriptPath = clustermapsScriptFile.canonicalPath
def spieceasiScriptFile = new File("${projectDir}/processes/spieceasi/run_spieceasi.R")
if( !spieceasiScriptFile.exists() ) {
    exit 1, "run_spieceasi.R not found in project directory"
}
def spieceasiScriptPath = spieceasiScriptFile.canonicalPath
def networkModulesScriptFile = new File("${projectDir}/processes/network_modules/network_modules.R")
if( !networkModulesScriptFile.exists() ) {
    exit 1, "network_modules.R not found in project directory"
}
def networkModulesScriptPath = networkModulesScriptFile.canonicalPath
def graphNetworkScriptFile = new File("${projectDir}/processes/graph_network/graph_network.py")
if( !graphNetworkScriptFile.exists() ) {
    exit 1, "graph_network.py not found in project directory"
}
def graphNetworkScriptPath = graphNetworkScriptFile.canonicalPath
def genomeCooccurrencePrepareScriptFile = new File("${projectDir}/processes/genome_cooccurrence/prepare_genome_cooccurrence.py")
def genomeCooccurrenceInferenceScriptFile = new File("${projectDir}/processes/genome_cooccurrence/infer_genome_proportionality.py")
def genomeCooccurrencePlotScriptFile = new File("${projectDir}/processes/genome_cooccurrence/plot_genome_cooccurrence.py")
for( scriptFile in [genomeCooccurrencePrepareScriptFile, genomeCooccurrenceInferenceScriptFile, genomeCooccurrencePlotScriptFile] ) {
    if( !scriptFile.exists() ) {
        exit 1, "Genome co-occurrence script not found: ${scriptFile}"
    }
}
def genomeCooccurrencePrepareScriptPath = genomeCooccurrencePrepareScriptFile.canonicalPath
def genomeCooccurrenceInferenceScriptPath = genomeCooccurrenceInferenceScriptFile.canonicalPath
def genomeCooccurrencePlotScriptPath = genomeCooccurrencePlotScriptFile.canonicalPath
def genomeCooccurrenceScriptHash = "${fileMd5(genomeCooccurrencePrepareScriptFile)}:${fileMd5(genomeCooccurrenceInferenceScriptFile)}:${fileMd5(genomeCooccurrencePlotScriptFile)}"
def masterSummaryScriptFile = new File("${projectDir}/processes/master_summary/build_master_asv_summary.py")
if( !masterSummaryScriptFile.exists() ) {
    exit 1, "summary/build_master_asv_summary.py not found in project directory"
}
def masterSummaryScriptPath = masterSummaryScriptFile.canonicalPath
def asvMagLinkScriptFile = new File("${projectDir}/processes/asv_mag_link/asv_mag_barrnap_linker.py")
if( !asvMagLinkScriptFile.exists() ) {
    exit 1, "asv_mag_barrnap_linker.py not found in project directory"
}
def asvMagLinkScriptPath = asvMagLinkScriptFile.canonicalPath
def plotAsvMagLinkScriptFile = new File("${projectDir}/processes/asv_mag_link/plot_asv_mag_link.py")
if( !plotAsvMagLinkScriptFile.exists() ) {
    exit 1, "plot_asv_mag_link.py not found in project directory"
}
def plotAsvMagLinkScriptPath = plotAsvMagLinkScriptFile.canonicalPath
def asvMagNetworkScriptFile = new File("${projectDir}/processes/asv_mag_network/asv_mag_network.py")
if( !asvMagNetworkScriptFile.exists() ) {
    exit 1, "asv_mag_network.py not found in project directory"
}
def asvMagNetworkScriptPath = asvMagNetworkScriptFile.canonicalPath
def asvMagGroupingDiagnosticsScriptFile = new File("${projectDir}/processes/asv_mag_network/mag_grouping_diagnostics.py")
if( !asvMagGroupingDiagnosticsScriptFile.exists() ) {
    exit 1, "mag_grouping_diagnostics.py not found in project directory"
}
def asvMagNetworkScriptHash = "${fileMd5(asvMagNetworkScriptFile)}:${fileMd5(asvMagGroupingDiagnosticsScriptFile)}"
def asvMagCurtainsScriptFile = new File("${projectDir}/processes/asv_mag_curtains/asv_mag_curtains.py")
if( !asvMagCurtainsScriptFile.exists() ) {
    exit 1, "ASV-MAG curtain script not found: ${asvMagCurtainsScriptFile}"
}
def asvMagCurtainsScriptPath = asvMagCurtainsScriptFile.canonicalPath
def asvMagCurtainsScriptHash = fileMd5(asvMagCurtainsScriptFile)
def groupGuildFunctionScriptFile = new File("${projectDir}/processes/group_guild_function/group_guild_function.py")
if( !groupGuildFunctionScriptFile.exists() ) {
    exit 1, "group_guild_function.py not found in project directory"
}
def groupGuildFunctionScriptPath = groupGuildFunctionScriptFile.canonicalPath
def groupGuildFunctionScriptHash = fileMd5(groupGuildFunctionScriptFile)
def moduleMagAnchorsScriptFile = new File("${projectDir}/processes/module_mag_anchors/summarize_module_mag_anchors.py")
if( !moduleMagAnchorsScriptFile.exists() ) {
    exit 1, "summarize_module_mag_anchors.py not found in project directory"
}
def moduleMagAnchorsScriptPath = moduleMagAnchorsScriptFile.canonicalPath
def powerAnalysisScriptFile = new File("${projectDir}/processes/power_analysis_pipeline/run_power_analysis_pipeline.sh")
if( !powerAnalysisScriptFile.exists() ) {
    exit 1, "run_power_analysis_pipeline.sh not found in project directory"
}
def powerAnalysisScriptPath = powerAnalysisScriptFile.canonicalPath
def brayPatientAwareScriptFile = new File("${projectDir}/processes/bray_patient_aware/run_bray_permanova_patient_aware.R")
if( !brayPatientAwareScriptFile.exists() ) {
    exit 1, "run_bray_permanova_patient_aware.R not found in project directory"
}
def brayPatientAwareScriptPath = brayPatientAwareScriptFile.canonicalPath
def plotBrayPatientAwareScriptFile = new File("${projectDir}/processes/bray_patient_aware/plot_bray_permanova_patient_aware.py")
if( !plotBrayPatientAwareScriptFile.exists() ) {
    exit 1, "plot_bray_permanova_patient_aware.py not found in project directory"
}
def plotBrayPatientAwareScriptPath = plotBrayPatientAwareScriptFile.canonicalPath
def taxonomicAbundanceObservedScriptFile = new File("${projectDir}/processes/taxonomy_patient_aware/run_taxonomic_abundance_analysis.py")
if( !taxonomicAbundanceObservedScriptFile.exists() ) {
    exit 1, "run_taxonomic_abundance_analysis.py not found in project directory"
}
def taxonomicAbundanceObservedScriptPath = taxonomicAbundanceObservedScriptFile.canonicalPath
def taxonomicSampleTypeObservedScriptFile = new File("${projectDir}/processes/taxonomy_patient_aware/run_taxonomic_sample_type_analysis.py")
if( !taxonomicSampleTypeObservedScriptFile.exists() ) {
    exit 1, "run_taxonomic_sample_type_analysis.py not found in project directory"
}
def taxonomicSampleTypeObservedScriptPath = taxonomicSampleTypeObservedScriptFile.canonicalPath
def plotTaxonomicObservedScriptFile = new File("${projectDir}/processes/taxonomy_patient_aware/plot_taxonomic_observed_analysis.py")
if( !plotTaxonomicObservedScriptFile.exists() ) {
    exit 1, "plot_taxonomic_observed_analysis.py not found in project directory"
}
def plotTaxonomicObservedScriptPath = plotTaxonomicObservedScriptFile.canonicalPath
def prepareLungStatusScriptFile = new File("${projectDir}/processes/lung_status_analysis/prepare_lung_status_data.py")
if( !prepareLungStatusScriptFile.exists() ) {
    exit 1, "prepare_lung_status_data.py not found in project directory"
}
def prepareLungStatusScriptPath = prepareLungStatusScriptFile.canonicalPath
def lungStatusAnalysisScriptFile = new File("${projectDir}/processes/lung_status_analysis/run_lung_status_analysis.R")
if( !lungStatusAnalysisScriptFile.exists() ) {
    exit 1, "run_lung_status_analysis.R not found in project directory"
}
def lungStatusAnalysisScriptPath = lungStatusAnalysisScriptFile.canonicalPath
def plotLungStatusScriptFile = new File("${projectDir}/processes/lung_status_analysis/plot_lung_status_analysis.py")
if( !plotLungStatusScriptFile.exists() ) {
    exit 1, "plot_lung_status_analysis.py not found in project directory"
}
def plotLungStatusScriptPath = plotLungStatusScriptFile.canonicalPath
def plotVocCorrScriptFile = new File("${projectDir}/processes/voc_correlation/plot_voc_corr.py")
if( !plotVocCorrScriptFile.exists() ) {
    exit 1, "plot_voc_corr.py not found in project directory"
}
def plotVocCorrScriptPath = plotVocCorrScriptFile.canonicalPath
def plotVocCorrScriptHash = fileMd5(plotVocCorrScriptFile)
def measurementAssociationScriptFile = new File("${projectDir}/processes/measurement_association/measurement_association.py")
if( !measurementAssociationScriptFile.exists() ) {
    exit 1, "measurement_association.py not found in project directory"
}
def measurementAssociationScriptPath = measurementAssociationScriptFile.canonicalPath
def measurementAssociationScriptHash = fileMd5(measurementAssociationScriptFile)
def moduleMeasurementAssociationScriptFile = new File("${projectDir}/processes/measurement_association/module_measurement_association.py")
if( !moduleMeasurementAssociationScriptFile.exists() ) {
    exit 1, "module_measurement_association.py not found in project directory"
}
def moduleMeasurementAssociationScriptPath = moduleMeasurementAssociationScriptFile.canonicalPath
def moduleMeasurementAssociationScriptHash = fileMd5(moduleMeasurementAssociationScriptFile)
def measurementAssociationRScriptFile = new File("${projectDir}/processes/measurement_association/run_measurement_association.R")
if( !measurementAssociationRScriptFile.exists() ) {
    exit 1, "run_measurement_association.R not found in project directory"
}
def measurementAssociationRScriptPath = measurementAssociationRScriptFile.canonicalPath
def measurementAssociationRScriptHash = fileMd5(measurementAssociationRScriptFile)
def titanPrepareScriptFile = new File("${projectDir}/processes/titan/prepare_titan.py")
def titanInstallScriptFile = new File("${projectDir}/processes/titan/install_titan2.R")
def titanAnalysisScriptFile = new File("${projectDir}/processes/titan/run_titan.R")
def titanCollectScriptFile = new File("${projectDir}/processes/titan/collect_titan.py")
def titanPlotScriptFile = new File("${projectDir}/processes/titan/plot_titan.py")
[titanPrepareScriptFile, titanInstallScriptFile, titanAnalysisScriptFile, titanCollectScriptFile, titanPlotScriptFile].each { scriptFile ->
    if( !scriptFile.exists() ) {
        exit 1, "TITAN script not found: ${scriptFile}"
    }
}
def titanPrepareScriptPath = titanPrepareScriptFile.canonicalPath
def titanInstallScriptPath = titanInstallScriptFile.canonicalPath
def titanAnalysisScriptPath = titanAnalysisScriptFile.canonicalPath
def titanCollectScriptPath = titanCollectScriptFile.canonicalPath
def titanPlotScriptPath = titanPlotScriptFile.canonicalPath
def titanPrepareScriptHash = fileMd5(titanPrepareScriptFile)
def titanInstallScriptHash = fileMd5(titanInstallScriptFile)
def titanAnalysisScriptHash = fileMd5(titanAnalysisScriptFile)
def titanCollectScriptHash = fileMd5(titanCollectScriptFile)
def titanPlotScriptHash = fileMd5(titanPlotScriptFile)

def microbialCompartmentPrepareScriptFile = new File("${projectDir}/processes/microbial_compartments/prepare_microbial_compartments.py")
def microbialCompartmentInferenceScriptFile = new File("${projectDir}/processes/microbial_compartments/infer_microbial_compartments.R")
def microbialCompartmentPosthocScriptFile = new File("${projectDir}/processes/microbial_compartments/posthoc_microbial_compartments.py")
def microbialCompartmentPlotScriptFile = new File("${projectDir}/processes/microbial_compartments/plot_microbial_compartments.py")
[microbialCompartmentPrepareScriptFile, microbialCompartmentInferenceScriptFile,
 microbialCompartmentPosthocScriptFile, microbialCompartmentPlotScriptFile].each { scriptFile ->
    if( !scriptFile.exists() ) exit 1, "Microbial-compartment script not found: ${scriptFile}"
}
def microbialCompartmentPrepareScriptPath = microbialCompartmentPrepareScriptFile.canonicalPath
def microbialCompartmentInferenceScriptPath = microbialCompartmentInferenceScriptFile.canonicalPath
def microbialCompartmentPosthocScriptPath = microbialCompartmentPosthocScriptFile.canonicalPath
def microbialCompartmentPlotScriptPath = microbialCompartmentPlotScriptFile.canonicalPath
def microbialCompartmentPrepareScriptHash = fileMd5(microbialCompartmentPrepareScriptFile)
def microbialCompartmentInferenceScriptHash = fileMd5(microbialCompartmentInferenceScriptFile)
def microbialCompartmentPosthocScriptHash = fileMd5(microbialCompartmentPosthocScriptFile)
def microbialCompartmentPlotScriptHash = fileMd5(microbialCompartmentPlotScriptFile)
def microbialStateInterpretationScriptFile = new File("${projectDir}/processes/microbial_state_interpretation/microbial_state_interpretation.py")
if( !microbialStateInterpretationScriptFile.exists() ) exit 1, "Microbial-state interpretation script not found: ${microbialStateInterpretationScriptFile}"
def microbialStateInterpretationScriptPath = microbialStateInterpretationScriptFile.canonicalPath
def microbialStateInterpretationScriptHash = fileMd5(microbialStateInterpretationScriptFile)
def ecologicalContextAtlasScriptFile = new File("${projectDir}/processes/ecological_context_atlas/ecological_context_atlas.py")
if( !ecologicalContextAtlasScriptFile.exists() ) exit 1, "Ecological-context atlas script not found: ${ecologicalContextAtlasScriptFile}"
def ecologicalContextAtlasScriptPath = ecologicalContextAtlasScriptFile.canonicalPath
def ecologicalContextAtlasScriptHash = "${fileMd5(ecologicalContextAtlasScriptFile)}:${fileMd5(asvMagNetworkScriptFile)}"
def communityTurnoverScriptFile = new File("${projectDir}/processes/community_turnover/community_turnover.py")
def communityTurnoverLcbdScriptFile = new File("${projectDir}/processes/community_turnover/lcbd_analysis.R")
def communityTurnoverPlotScriptFile = new File("${projectDir}/processes/community_turnover/plot_community_turnover.py")
[communityTurnoverScriptFile, communityTurnoverLcbdScriptFile, communityTurnoverPlotScriptFile].each { scriptFile ->
    if( !scriptFile.exists() ) exit 1, "Community-turnover script not found: ${scriptFile}"
}
def communityTurnoverScriptPath = communityTurnoverScriptFile.canonicalPath
def communityTurnoverLcbdScriptPath = communityTurnoverLcbdScriptFile.canonicalPath
def communityTurnoverPlotScriptPath = communityTurnoverPlotScriptFile.canonicalPath
def communityTurnoverScriptHash = fileMd5(communityTurnoverScriptFile)
def communityTurnoverLcbdScriptHash = fileMd5(communityTurnoverLcbdScriptFile)
def communityTurnoverPlotScriptHash = fileMd5(communityTurnoverPlotScriptFile)
def groupingDiagnosticsScriptFile = new File("${projectDir}/processes/grouping_diagnostics/grouping_diagnostics.py")
if( !groupingDiagnosticsScriptFile.exists() ) {
    exit 1, "grouping_diagnostics.py not found in project directory"
}
def groupingDiagnosticsScriptPath = groupingDiagnosticsScriptFile.canonicalPath
def groupingDiagnosticsScriptHash = fileMd5(groupingDiagnosticsScriptFile)
def communityPredictorScriptFile = new File("${projectDir}/processes/community_predictor_comparison/community_predictor_comparison.py")
if( !communityPredictorScriptFile.exists() ) {
    exit 1, "community_predictor_comparison.py not found in project directory"
}
def communityPredictorScriptPath = communityPredictorScriptFile.canonicalPath
def communityPredictorScriptHash = fileMd5(communityPredictorScriptFile)
def groupLabelAugmentationScriptFile = new File("${projectDir}/processes/group_label_augmentation/group_label_augmentation.py")
if( !groupLabelAugmentationScriptFile.exists() ) {
    exit 1, "group_label_augmentation.py not found in project directory"
}
def groupLabelAugmentationScriptPath = groupLabelAugmentationScriptFile.canonicalPath
def groupLabelAugmentationScriptHash = fileMd5(groupLabelAugmentationScriptFile)
def emptyModulesScriptFile = new File("${projectDir}/processes/master_summary/empty_modules.tsv")
if( !emptyModulesScriptFile.exists() ) {
    exit 1, "empty_modules.tsv not found in project directory"
}
def emptyModulesPath = emptyModulesScriptFile.canonicalPath
def sinaConfig = config.sina ?: [:]
def sinaDownloadSubdir = (sinaConfig.download_subdir ?: 'sina_reference').toString()
def sinaDownloadDir = new File(outputDir, sinaDownloadSubdir)
sinaDownloadDir.mkdirs()
def defaultSinaReferenceFile = new File(sinaDownloadDir, sinaReferenceFilename)
def configuredSinaReference = sinaConfig.reference ? resolveOptionalPath(sinaConfig.reference, configRoot) : null
def initialSinaReferencePath = configuredSinaReference ?: defaultSinaReferenceFile.canonicalPath
File sinaReferenceFile = new File(initialSinaReferencePath)
if( aspirePhase != 'analysis' && !sinaReferenceFile.exists() ) {
    if( configuredSinaReference && sinaReferenceFile.canonicalPath != defaultSinaReferenceFile.canonicalPath ) {
        log.warn "Configured SINA reference not found at ${sinaReferenceFile.canonicalPath}; downloading to ${defaultSinaReferenceFile.canonicalPath}."
        sinaReferenceFile = defaultSinaReferenceFile
    }
    def referenceUrl = sinaConfig.reference_url ?: defaultSinaReferenceUrl
    log.info "Downloading SILVA reference (${sinaReferenceFilename}) from ${referenceUrl}"
    downloadReference(referenceUrl, sinaReferenceFile)
}
if( aspirePhase != 'analysis' && !sinaReferenceFile.exists() ) {
    exit 1, "Failed to obtain SINA reference file at ${sinaReferenceFile}"
}
def sinaReferencePath = sinaReferenceFile.canonicalPath
def sinaRegionsRaw = sinaConfig.regions
List<String> sinaRegionList
if( sinaRegionsRaw instanceof List ) {
    sinaRegionList = sinaRegionsRaw.collect { it.toString().trim() }.findAll { it }
} else if( sinaRegionsRaw ) {
    sinaRegionList = sinaRegionsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
} else {
    sinaRegionList = ['V4']
}
if( !sinaRegionList ) {
    sinaRegionList = ['V4']
}
def sinaRegionsArg = sinaRegionList.join(',')
def sinaTrimTarget = sinaConfig.trim_to ?: sinaRegionList[0]
def sinaBatchSize = sinaConfig.batch_size ? (sinaConfig.batch_size as int) : 1000000
def sinaKeepGaps = (sinaConfig.keep_gaps ?: false) as boolean
def sinaVerbose = (sinaConfig.verbose ?: false) as boolean
def sinaIdColumn = sinaConfig.id_column ?: 'ASV_ID'
def sinaThreads = sinaConfig.threads ? (sinaConfig.threads as int) : pipelineThreads
def generalStatsConfig = config.general_stats ?: [:]
def generalStatsEnabled = generalStatsConfig.containsKey('enabled') ? (generalStatsConfig.enabled as boolean) : true
def rawInputFiles = generalStatsEnabled ? sampleRecords.collectMany { rec ->
    def list = [rec.r1]
    if( rec.paired && rec.r2 ) {
        list << rec.r2
    }
    return list
}.unique() : []
def fastpOutputFiles = generalStatsEnabled ? sampleRecords.collectMany { rec ->
    def outputs = [
        new File(dirMap.fastp, "${rec.sample_id}_R1.fastq.gz").absolutePath
    ]
    outputs << new File(dirMap.fastp, "${rec.sample_id}_R2.fastq.gz").absolutePath
    return outputs
} : []
def filteredOutputFiles = generalStatsEnabled ? sampleRecords.collect { rec ->
    new File(dirMap.filter, "${rec.sample_id}.filtered.fasta.gz").absolutePath
} : []
def generalStatsRawArgs = joinShellArgs(rawInputFiles)
def generalStatsFastpArgs = joinShellArgs(fastpOutputFiles)
def generalStatsFilteredArgs = joinShellArgs(filteredOutputFiles)
def taxonomyScriptFile = new File("${projectDir}/processes/taxonomy/qiime_vs_classifier.py")
if( !taxonomyScriptFile.exists() ) {
    exit 1, "qiime_vs_classifier.py not found in project directory"
}
def taxonomyScriptPath = taxonomyScriptFile.canonicalPath
def taxonomyConfig = config.taxonomy ?: [:]
def taxonomyDownloadSubdir = (taxonomyConfig.download_subdir ?: 'taxonomy_reference').toString()
def taxonomyDownloadDir = new File(outputDir, taxonomyDownloadSubdir)
taxonomyDownloadDir.mkdirs()
def taxonomyTaxFilename = taxonomyConfig.ref_taxonomy_filename ?: defaultTaxonomyTaxFilename
def taxonomySeqFilename = taxonomyConfig.ref_sequences_filename ?: defaultTaxonomySeqsFilename

def taxonomyRefTaxonomyResolved = taxonomyConfig.ref_taxonomy ? resolveOptionalPath(taxonomyConfig.ref_taxonomy, configRoot) : null
def taxonomyRefTaxonomyFile = taxonomyRefTaxonomyResolved ? new File(taxonomyRefTaxonomyResolved) : new File(taxonomyDownloadDir, taxonomyTaxFilename)
if( aspirePhase != 'analysis' && !taxonomyRefTaxonomyFile.exists() ) {
    def taxonomyUrl = taxonomyConfig.ref_taxonomy_url ?: defaultTaxonomyTaxUrl
    if( !taxonomyUrl ) {
        exit 1, "taxonomy.ref_taxonomy or taxonomy.ref_taxonomy_url must be provided to obtain SILVA taxonomy artifact"
    }
    log.info "Downloading SILVA taxonomy artifact from ${taxonomyUrl}"
    downloadReference(taxonomyUrl, taxonomyRefTaxonomyFile)
}
if( aspirePhase != 'analysis' && !taxonomyRefTaxonomyFile.exists() ) {
    exit 1, "Taxonomy reference taxonomy file not found: ${taxonomyRefTaxonomyFile}"
}

def taxonomyRefSequencesResolved = taxonomyConfig.ref_sequences ? resolveOptionalPath(taxonomyConfig.ref_sequences, configRoot) : null
def taxonomyRefSequencesFile = taxonomyRefSequencesResolved ? new File(taxonomyRefSequencesResolved) : new File(taxonomyDownloadDir, taxonomySeqFilename)
if( aspirePhase != 'analysis' && !taxonomyRefSequencesFile.exists() ) {
    def taxonomySeqsUrl = taxonomyConfig.ref_sequences_url ?: defaultTaxonomySeqsUrl
    if( !taxonomySeqsUrl ) {
        exit 1, "taxonomy.ref_sequences or taxonomy.ref_sequences_url must be provided to obtain SILVA sequences artifact"
    }
    log.info "Downloading SILVA sequences artifact from ${taxonomySeqsUrl}"
    downloadReference(taxonomySeqsUrl, taxonomyRefSequencesFile)
}
if( aspirePhase != 'analysis' && !taxonomyRefSequencesFile.exists() ) {
    exit 1, "Taxonomy reference sequences file not found: ${taxonomyRefSequencesFile}"
}

def taxonomyRefTaxonomy = taxonomyRefTaxonomyFile.canonicalPath
def taxonomyRefSequences = taxonomyRefSequencesFile.canonicalPath
def taxonomyOutputName = taxonomyConfig.output_tsv ?: 'ASV_SILVA_tax.full-length.vsearch.tsv'
def taxonomyStatsName = taxonomyConfig.stats_tsv ?: 'ASV_SILVA_stats.full-length.vsearch.tsv'
def taxonomyUppercaseName = taxonomyConfig.uppercase_fasta ?: 'ASVs.upper.fasta'
def taxonomyUppercasePlainName = taxonomyUppercaseName.toString().endsWith('.gz') ? taxonomyUppercaseName.toString()[0..-4] : taxonomyUppercaseName.toString()
def taxonomyUppercaseGzName = taxonomyUppercaseName.toString().endsWith('.gz') ? taxonomyUppercaseName.toString() : "${taxonomyUppercasePlainName}.gz"
def taxonomyThreads = taxonomyConfig.threads ? (taxonomyConfig.threads as int) : pipelineThreads
def mitoConfig = config.mito ?: [:]
def mitoEnabled = mitoConfig.containsKey('enabled') ? (mitoConfig.enabled as boolean) : true
def defaultMitoChunkDir = new File(dirMap.asv, "chunks").canonicalPath
def mitoChunkDirPath = mitoConfig.chunk_dir ? resolveOutputRelative(mitoConfig.chunk_dir.toString(), outputDir) : defaultMitoChunkDir
def defaultMitoOutputDir = new File(dirMap.mito, "mitomap").canonicalPath
def mitoOutputDirPath = mitoConfig.output_dir ? resolveOutputRelative(mitoConfig.output_dir.toString(), outputDir) : defaultMitoOutputDir
def mitoBlastDbPath = resolveOptionalPath(mitoConfig.mito_db ?: 'ref_db/mito_ncbi', configRoot)
def mitoBiofDbPath = resolveOptionalPath(mitoConfig.biof_db ?: 'ref_db/ssu_pipeline_contaminants', configRoot)
def mitoBlastFastaPath = mitoConfig.mito_fasta ? resolveOptionalPath(mitoConfig.mito_fasta.toString(), configRoot) : null
def mitoBiofFastaPath = mitoConfig.contaminant_fasta ? resolveOptionalPath(mitoConfig.contaminant_fasta.toString(), configRoot) :
    (mitoConfig.biof_fasta ? resolveOptionalPath(mitoConfig.biof_fasta.toString(), configRoot) : null)
def mitoChunkSize = mitoConfig.chunk_size ? (mitoConfig.chunk_size as int) : 10
def mitomasterWorkers = mitoConfig.mitomaster_workers ? (mitoConfig.mitomaster_workers as int) : 8
def mitomasterRetries = mitoConfig.mitomaster_retries ? (mitoConfig.mitomaster_retries as int) : 4
def mitomasterTimeout = mitoConfig.mitomaster_timeout ? (mitoConfig.mitomaster_timeout as int) : 90
def mitomasterHeaderMode = mitoConfig.mitomaster_header_mode ?: 'first'
def mitoRunMitomaster = mitoConfig.containsKey('run_mitomaster') ? (mitoConfig.run_mitomaster as boolean) : true
def mitoBlastThreads = mitoConfig.blast_threads ? (mitoConfig.blast_threads as int) : pipelineThreads
def mitoPrefix = mitoConfig.prefix ?: 'nontarget'
def mitoFormats = mitoConfig.formats ?: 'svg,pdf'
def mitoMinPident = mitoConfig.min_pident != null ? (mitoConfig.min_pident as double) : 97.0
def mitoMinPercov = mitoConfig.min_percov != null ? (mitoConfig.min_percov as double) : 51.0
def mitoMitoSubstring = mitoConfig.mitochondria_substring ?: 'mitochondria'
def mitoFeatureCol = mitoConfig.feature_col ?: 'Feature ID'
def mitoTaxonCol = mitoConfig.taxon_col ?: 'Taxon'
def mitoConsensusCol = mitoConfig.consensus_col ?: 'Consensus'
def mitoSteps = mitoConfig.steps ?: 'BioFactorial,Qiime_NB_FULL,MITOMASTER,BLAST_mito'
def mitoHostFirstStep = mitoConfig.host_first_step ?: 'BioFactorial'
def mitoFigsize = mitoConfig.figsize ?: '10x6'
def mitoStyle = mitoConfig.style ?: 'whitegrid'
def mitoDpi = mitoConfig.dpi ? (mitoConfig.dpi as int) : 300
def mitoNoPlots = (mitoConfig.no_plots ?: false) as boolean
def filterCountsConfig = config.filter_counts ?: [:]
def filterCountsEnabled = filterCountsConfig.containsKey('enabled') ? (filterCountsConfig.enabled as boolean) : true
if( filterCountsEnabled && !mitoEnabled ) {
    exit 1, "filter_counts.enabled requires mito.enabled to be true"
}
def filterCountsMetadataPath = filterCountsConfig.metadata ? resolveOptionalPath(filterCountsConfig.metadata, configRoot) : null
if( filterCountsEnabled && filterCountsMetadataPath && !new File(filterCountsMetadataPath).exists() ) {
    exit 1, "filter_counts.metadata not found: ${filterCountsMetadataPath}"
}
def filterCountsOutputName = filterCountsConfig.output ?: 'ASV_target.tsv'
def filterCountsGroupCol = filterCountsConfig.group_col ?: 'Depth'
def filterCountsMinGroup = filterCountsConfig.min_group_size ? (filterCountsConfig.min_group_size as int) : 3
def filterCountsAbundance = filterCountsConfig.abundance_threshold != null ? (filterCountsConfig.abundance_threshold as double) : 0.005d
def filterCountsSampleCol = filterCountsConfig.sample_id_col ?: 'longID'
def filterCountsMinConsensus = filterCountsConfig.min_consensus != null ? (filterCountsConfig.min_consensus as double) : 0d
def filterCountsTaxonCol = filterCountsConfig.taxon_col ?: 'Taxon'
def filterCountsConsensusCol = filterCountsConfig.consensus_col ?: 'Consensus'
def filterCountsBiofactorialCol = filterCountsConfig.biofactorial_col ?: 'BioFactorial'
def filterCountsMitoColsRaw = filterCountsConfig.mito_cols
List<String> filterCountsMitoCols
if( filterCountsMitoColsRaw instanceof List ) {
    filterCountsMitoCols = filterCountsMitoColsRaw.collect { it.toString() }
} else if( filterCountsMitoColsRaw ) {
    filterCountsMitoCols = filterCountsMitoColsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
} else {
    filterCountsMitoCols = ['MITOMASTER','BLAST_mito']
}
def filterCountsExcludeTaxaRaw = filterCountsConfig.exclude_taxa
List<String> filterCountsExcludeTaxa = []
if( filterCountsExcludeTaxaRaw instanceof List ) {
    filterCountsExcludeTaxa = filterCountsExcludeTaxaRaw.collect { item ->
        if( item instanceof Map ) {
            def rank = item.rank ?: item.level
            def value = item.value ?: item.taxon ?: item.name
            (rank && value) ? "${rank}:${value}".toString() : ''
        } else {
            item.toString()
        }
    }.findAll { it?.trim() }
} else if( filterCountsExcludeTaxaRaw instanceof Map ) {
    filterCountsExcludeTaxaRaw.each { rank, values ->
        if( values instanceof List ) {
            values.each { value -> filterCountsExcludeTaxa << "${rank}:${value}".toString() }
        } else if( values ) {
            filterCountsExcludeTaxa << "${rank}:${values}".toString()
        }
    }
} else if( filterCountsExcludeTaxaRaw ) {
    filterCountsExcludeTaxa = filterCountsExcludeTaxaRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
}
def filterCountsSaveIntermediates = (filterCountsConfig.save_intermediates ?: false) as boolean
def defaultFilterMitoDir = new File(dirMap.mito, "ASVs").canonicalPath
def filterCountsMitoDir = filterCountsConfig.mito_output_dir ? resolveOutputRelative(filterCountsConfig.mito_output_dir.toString(), outputDir) : defaultFilterMitoDir
if( mitoEnabled ) {
    ensureBlastReferenceExists(mitoBlastDbPath, mitoBlastFastaPath, 'mitochondrial')
    ensureBlastReferenceExists(mitoBiofDbPath, mitoBiofFastaPath, 'contaminant')
    new File(mitoChunkDirPath).parentFile?.mkdirs()
    new File(mitoOutputDirPath).mkdirs()
}
def sankeyConfig = config.sankey ?: [:]
def sankeyEnabled = (sankeyConfig.enabled ?: false) as boolean
def sankeyMetadataPath = sankeyConfig.metadata ? resolveOptionalPath(sankeyConfig.metadata, configRoot) : resolveOptionalPath("ref_db/asv_cruise_metadata.tsv", configRoot)
if( sankeyEnabled && (!sankeyMetadataPath || !new File(sankeyMetadataPath).exists()) ) {
    exit 1, "Sankey metadata file not found: ${sankeyMetadataPath}"
}
def sankeySubDir = sankeyConfig.sub_dir ?: '.'
def sankeySampCol = sankeyConfig.sample_col ?: 'lmp_id'
def sankeyGroupCol = sankeyConfig.group1_col ?: 'group1'
def sankeyColorCol = sankeyConfig.color_col ?: 'Color'
if( sankeyEnabled ) {
    def sankeyPalettePath = sankeyConfig.palette_file ? resolveOptionalPath(sankeyConfig.palette_file.toString(), configRoot) : null
    def sankeyAssets = prepareMetadataAssets(
        sankeyMetadataPath,
        sankeySampCol.toString(),
        sankeyGroupCol.toString(),
        sankeyColorCol.toString(),
        sankeyPalettePath,
        new File(dirMap.metadata, "sankey_${safeFilename(sankeyGroupCol.toString())}")
    )
    sankeyMetadataPath = sankeyAssets.metadata
    log.info "Sankey metadata palette: ${sankeyAssets.palette}"
}
def sankeyKeepTypesRaw = sankeyConfig.keep_types
List<String> sankeyKeepTypes = []
if( sankeyKeepTypesRaw instanceof List ) {
    sankeyKeepTypes = sankeyKeepTypesRaw.collect { it.toString().trim() }.findAll { it }
} else if( sankeyKeepTypesRaw ) {
    sankeyKeepTypes = sankeyKeepTypesRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
}
def sankeyVerticalOrderRaw = sankeyConfig.vertical_order ?: sankeyConfig.group_order
List<String> sankeyVerticalOrder = normalizePresetList(sankeyVerticalOrderRaw, [], config.order_presets ?: [:])
def sankeyOutputPrefix = sankeyConfig.output_prefix ?: "metadata/data_loss_sankey"
def sankeyTitle = sankeyConfig.title ?: "Data Loss Flow"
def sankeyMakeLabeled = (sankeyConfig.make_labeled == null) ? true : (sankeyConfig.make_labeled as boolean)
def sankeyMakeUnlabeled = (sankeyConfig.make_unlabeled == null) ? true : (sankeyConfig.make_unlabeled as boolean)
def sankeyArrangement = (sankeyConfig.arrangement ?: 'snap').toString().trim().toLowerCase()
if( !(sankeyArrangement in ['snap', 'perpendicular', 'freeform', 'fixed']) ) {
    exit 1, "Invalid sankey.arrangement '${sankeyArrangement}'. Allowed: snap, perpendicular, freeform, fixed"
}
if( sankeyEnabled && filterCountsEnabled && !filterCountsSaveIntermediates ) {
    exit 1, "Sankey requires filter_counts.save_intermediates to be true to access intermediate tables"
}
if( sankeyEnabled && !filterCountsEnabled ) {
    exit 1, "Sankey requires filter_counts.enabled to be true"
}
if( sankeyEnabled && !generalStatsEnabled ) {
    exit 1, "Sankey requires general_stats.enabled to be true"
}
def analysisCfg = parseMetadataAndBasicAnalysisConfig(config, configRoot, outputDir, filterCountsEnabled as boolean, generalStatsEnabled as boolean, pipelineThreads as int, dirMap)
analysisCfg.each { key, value -> binding.setVariable(key as String, value) }

def asvTimeDepthCurtainConfig = config.asv_time_depth_curtain ?: [:]
boolean asvTimeDepthCurtainEnabled = (asvTimeDepthCurtainConfig.enabled ?: false) as boolean
def asvTimeDepthCurtainOutputDirAbs = new File(
    outputDir, (asvTimeDepthCurtainConfig.output_dir ?: 'asv_time_depth_curtain').toString()
).canonicalPath
def asvTimeDepthCurtainSampleCol = (asvTimeDepthCurtainConfig.sample_col ?: 'sampleID').toString()
def asvTimeDepthCurtainCruiseCol = (asvTimeDepthCurtainConfig.cruise_col ?: 'Cruise').toString()
def asvTimeDepthCurtainDateCol = (asvTimeDepthCurtainConfig.date_col ?: 'Date').toString()
def asvTimeDepthCurtainDepthCol = (asvTimeDepthCurtainConfig.depth_col ?: 'Depth').toString()
def asvTimeDepthCurtainGroupCol = (asvTimeDepthCurtainConfig.group_col ?: 'o2_subcompartment_final').toString()
def asvTimeDepthCurtainPalette = (asvTimeDepthCurtainConfig.palette ?: '').toString()
def asvTimeDepthCurtainOrderRaw = asvTimeDepthCurtainConfig.group_order ?: []
def asvTimeDepthCurtainOrder = asvTimeDepthCurtainOrderRaw instanceof List ?
    asvTimeDepthCurtainOrderRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
    asvTimeDepthCurtainOrderRaw.toString()
def asvTimeDepthCurtainExcludePattern = (asvTimeDepthCurtainConfig.exclude_label_pattern ?: '(?i)(?:other|outlier)').toString()
double asvTimeDepthCurtainMaximumDepth = (asvTimeDepthCurtainConfig.maximum_depth ?: 210.0) as double
double asvTimeDepthCurtainDepthStep = (asvTimeDepthCurtainConfig.depth_step_m ?: 1.0) as double
int asvTimeDepthCurtainTimeSubdivisions = (asvTimeDepthCurtainConfig.time_subdivisions_per_month ?: asvTimeDepthCurtainConfig.time_subdivisions ?: 4) as int
double asvTimeDepthCurtainTimeSigma = (asvTimeDepthCurtainConfig.time_sigma_months ?: asvTimeDepthCurtainConfig.time_sigma_cruises ?: 0.75) as double
double asvTimeDepthCurtainContourDepthSigma = (asvTimeDepthCurtainConfig.contour_visual_depth_sigma_m ?: 5.0) as double
def asvTimeDepthCurtainSourceMode = (asvTimeDepthCurtainConfig.source_mode ?: 'asv').toString().trim().toLowerCase()
def asvTimeDepthCurtainBasinGrid = asvTimeDepthCurtainConfig.basin_grid ?
    resolveOptionalPath(asvTimeDepthCurtainConfig.basin_grid.toString(), configRoot) : null
def asvTimeDepthCurtainBasinCells = asvTimeDepthCurtainConfig.basin_cells ?
    resolveOptionalPath(asvTimeDepthCurtainConfig.basin_cells.toString(), configRoot) : null
def asvTimeDepthCurtainRenewalEvents = asvTimeDepthCurtainConfig.renewal_events ?
    resolveOptionalPath(asvTimeDepthCurtainConfig.renewal_events.toString(), configRoot) : null
def asvTimeDepthCurtainRenewalDateCol = (asvTimeDepthCurtainConfig.renewal_date_col ?: 'start_date').toString()
double asvTimeDepthCurtainBasePointSize = (asvTimeDepthCurtainConfig.base_point_size ?: 7.5) as double
double asvTimeDepthCurtainAsvPointSize = (asvTimeDepthCurtainConfig.asv_point_size ?: 42.0) as double
double asvTimeDepthCurtainAsvPointEdgeWidth = (asvTimeDepthCurtainConfig.asv_point_edge_width ?: 1.2) as double
def asvTimeDepthCurtainFormatsRaw = asvTimeDepthCurtainConfig.formats ?: ['pdf','png','svg']
def asvTimeDepthCurtainFormats = asvTimeDepthCurtainFormatsRaw instanceof List ?
    asvTimeDepthCurtainFormatsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
    asvTimeDepthCurtainFormatsRaw.toString()
if( asvTimeDepthCurtainEnabled && !metadataPlotsEnabled ) {
    exit 1, "asv_time_depth_curtain.enabled requires metadata_plots.enabled"
}
if( asvTimeDepthCurtainEnabled && !asvTimeDepthCurtainPalette ) {
    exit 1, "asv_time_depth_curtain.palette must define label=color entries"
}
if( asvTimeDepthCurtainEnabled && !(asvTimeDepthCurtainSourceMode in ['asv', 'basin']) ) {
    exit 1, "asv_time_depth_curtain.source_mode must be 'asv' or 'basin'"
}
if( asvTimeDepthCurtainEnabled && asvTimeDepthCurtainSourceMode == 'basin' &&
    (!asvTimeDepthCurtainBasinGrid || !asvTimeDepthCurtainBasinCells) ) {
    exit 1, "BASIN curtain mode requires basin_grid and basin_cells"
}

def indicatorCfg = parseIndicatorAndNetworkConfig(config, configRoot, outputDir, pipelineThreads as int, analysisCfg)
indicatorCfg.each { key, value -> binding.setVariable(key as String, value) }

def parseMetadataAndBasicAnalysisConfig(config, File configRoot, String outputDir, boolean filterCountsEnabled, boolean generalStatsEnabled, int pipelineThreads, Map dirMap) {
    def metadataPlotsConfig = config.metadata_plots ?: [:]
    boolean metadataPlotsEnabled = metadataPlotsConfig.containsKey('enabled') ? (metadataPlotsConfig.enabled as boolean) : true
    if( metadataPlotsEnabled && !filterCountsEnabled ) {
        exit 1, "metadata_plots.enabled requires filter_counts.enabled to be true"
    }
    if( metadataPlotsEnabled && !generalStatsEnabled ) {
        exit 1, "metadata_plots.enabled requires general_stats.enabled to be true"
    }
    def metadataPlotsMetadataPath = metadataPlotsConfig.metadata ? resolveOptionalPath(metadataPlotsConfig.metadata, configRoot) : resolveOptionalPath("ref_db/asv_cruise_metadata.tsv", configRoot)
    if( metadataPlotsEnabled && (!metadataPlotsMetadataPath || !new File(metadataPlotsMetadataPath).exists()) ) {
        exit 1, "metadata_plots metadata file not found: ${metadataPlotsMetadataPath}"
    }
    def metadataPlotsSubDir = metadataPlotsConfig.sub_dir ?: '.'
    def metadataPlotsSampleCol = metadataPlotsConfig.sample_col ?: 'sampleID'
    def metadataPlotsTypeCol = metadataPlotsConfig.type_col ?: (metadataPlotsConfig.group1_col ?: 'Depth')
    def metadataPlotsColorCol = metadataPlotsConfig.color_col ?: 'Color'
    if( metadataPlotsEnabled ) {
        def configuredPalettePath = metadataPlotsConfig.palette_file ? resolveOptionalPath(metadataPlotsConfig.palette_file.toString(), configRoot) : null
        def metadataAssets = prepareMetadataAssets(
            metadataPlotsMetadataPath,
            metadataPlotsSampleCol.toString(),
            metadataPlotsTypeCol.toString(),
            metadataPlotsColorCol.toString(),
            configuredPalettePath,
            new File(dirMap.metadata, safeFilename(metadataPlotsTypeCol.toString()))
        )
        metadataPlotsMetadataPath = metadataAssets.metadata
        log.info "Metadata palette: ${metadataAssets.palette}"
    }
    def metadataPlotsBiochemAssignmentsPath = metadataPlotsConfig.biochem_assignments ? resolveOptionalPath(metadataPlotsConfig.biochem_assignments, configRoot) : null
    def metadataPlotsBiochemSampleCol = metadataPlotsConfig.biochem_sample_col ?: 'cruise_year_month_depth'
    def metadataPlotsStratificationTimeseriesPath = metadataPlotsConfig.stratification_timeseries ? resolveOptionalPath(metadataPlotsConfig.stratification_timeseries, configRoot) : null
    def metadataPlotsStratMetaJoinCol = metadataPlotsConfig.strat_meta_join_col ?: 'Cruise'
    def metadataPlotsStratJoinCol = metadataPlotsConfig.strat_join_col ?: 'Cruise'
    def metadataPlotsCruiseGroupPath = metadataPlotsConfig.cruise_group_assignments ? resolveOptionalPath(metadataPlotsConfig.cruise_group_assignments, configRoot) : null
    def metadataPlotsCruiseGroupMetaJoinCol = metadataPlotsConfig.cruise_group_meta_join_col ?: 'Cruise'
    def metadataPlotsCruiseGroupJoinCol = metadataPlotsConfig.cruise_group_join_col ?: 'Cruise'
    def metadataCruiseGroupIncludeRaw = metadataPlotsConfig.cruise_group_include_cols ?: ['cruise_group','max_prob','resp_entropy_normalized','assignment_uncertain','PC1','PC2']
    List<String> metadataPlotsCruiseGroupIncludeCols = metadataCruiseGroupIncludeRaw instanceof List ?
        metadataCruiseGroupIncludeRaw.collect { it.toString().trim() }.findAll { it } :
        metadataCruiseGroupIncludeRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    def metadataGroupNormalizationConfig = metadataPlotsConfig.group_normalization instanceof Map ? metadataPlotsConfig.group_normalization : [:]
    boolean metadataGroupNormalizationEnabled = metadataGroupNormalizationConfig.containsKey('enabled') ? (metadataGroupNormalizationConfig.enabled as boolean) : false
    def metadataGroupNormalizationColsRaw = metadataGroupNormalizationConfig.columns ?: []
    List<String> metadataGroupNormalizationCols = metadataGroupNormalizationColsRaw instanceof List ?
        metadataGroupNormalizationColsRaw.collect { it.toString().trim() }.findAll { it } :
        metadataGroupNormalizationColsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    def metadataGroupNormalizationPattern = metadataGroupNormalizationConfig.pattern ? metadataGroupNormalizationConfig.pattern.toString() : ''
    def metadataGroupNormalizationReplacement = metadataGroupNormalizationConfig.replacement ? metadataGroupNormalizationConfig.replacement.toString().trim() : 'outlier'
    boolean metadataGroupNormalizationPreserveSource = metadataGroupNormalizationConfig.containsKey('preserve_source') ? (metadataGroupNormalizationConfig.preserve_source as boolean) : true
    if( metadataGroupNormalizationEnabled && (metadataGroupNormalizationCols.isEmpty() || !metadataGroupNormalizationPattern) ) {
        exit 1, "metadata_plots.group_normalization requires columns and pattern when enabled"
    }
    def metadataBiochemIncludeRaw = metadataPlotsConfig.biochem_include_cols
    List<String> metadataPlotsBiochemIncludeCols = []
    if( metadataBiochemIncludeRaw instanceof List ) {
        metadataPlotsBiochemIncludeCols = metadataBiochemIncludeRaw.collect { it.toString().trim() }.findAll { it }
    } else if( metadataBiochemIncludeRaw ) {
        metadataPlotsBiochemIncludeCols = metadataBiochemIncludeRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def metadataStratIncludeRaw = metadataPlotsConfig.strat_include_cols
    List<String> metadataPlotsStratIncludeCols = []
    if( metadataStratIncludeRaw instanceof List ) {
        metadataPlotsStratIncludeCols = metadataStratIncludeRaw.collect { it.toString().trim() }.findAll { it }
    } else if( metadataStratIncludeRaw ) {
        metadataPlotsStratIncludeCols = metadataStratIncludeRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def metadataBiochemMetaJoinRaw = metadataPlotsConfig.biochem_meta_join_cols
    List<String> metadataPlotsBiochemMetaJoinCols = []
    if( metadataBiochemMetaJoinRaw instanceof List ) {
        metadataPlotsBiochemMetaJoinCols = metadataBiochemMetaJoinRaw.collect { it.toString().trim() }.findAll { it }
    } else if( metadataBiochemMetaJoinRaw ) {
        metadataPlotsBiochemMetaJoinCols = metadataBiochemMetaJoinRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def metadataBiochemJoinRaw = metadataPlotsConfig.biochem_join_cols
    List<String> metadataPlotsBiochemJoinCols = []
    if( metadataBiochemJoinRaw instanceof List ) {
        metadataPlotsBiochemJoinCols = metadataBiochemJoinRaw.collect { it.toString().trim() }.findAll { it }
    } else if( metadataBiochemJoinRaw ) {
        metadataPlotsBiochemJoinCols = metadataBiochemJoinRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    if( metadataPlotsBiochemMetaJoinCols.size() != metadataPlotsBiochemJoinCols.size() ) {
        exit 1, "metadata_plots.biochem_meta_join_cols and metadata_plots.biochem_join_cols must have the same number of entries"
    }
    def metadataKeepTypesRaw = metadataPlotsConfig.keep_types
    List<String> metadataKeepTypes = []
    if( metadataKeepTypesRaw instanceof List ) {
        metadataKeepTypes = metadataKeepTypesRaw.collect { it.toString().trim() }.findAll { it }
    } else if( metadataKeepTypesRaw ) {
        metadataKeepTypes = metadataKeepTypesRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def metadataPlotsSubtractionCol = metadataPlotsConfig.subtraction_group_col ?: metadataPlotsTypeCol
    def metadataSubtractionGroupsRaw = metadataPlotsConfig.subtraction_groups
    List<String> metadataPlotsSubtractionGroups = []
    if( metadataSubtractionGroupsRaw instanceof List ) {
        metadataPlotsSubtractionGroups = metadataSubtractionGroupsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( metadataSubtractionGroupsRaw ) {
        metadataPlotsSubtractionGroups = metadataSubtractionGroupsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def metadataGroupOrderRaw = metadataPlotsConfig.group_order
    List<String> metadataPlotsGroupOrder = normalizePresetList(metadataGroupOrderRaw, [], config.order_presets ?: [:])
    def metadataIncludeRankRaw = metadataPlotsConfig.include_rank
    List<String> metadataIncludeRank = []
    if( metadataIncludeRankRaw instanceof List ) {
        metadataIncludeRank = metadataIncludeRankRaw.collect { it.toString().trim() }.findAll { it }
    } else if( metadataIncludeRankRaw ) {
        metadataIncludeRank = metadataIncludeRankRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def metadataPlotsMitoThreshold = metadataPlotsConfig.mito_threshold_line != null ? (metadataPlotsConfig.mito_threshold_line as double) : 1000d
    boolean metadataPlotsRunMicro = metadataPlotsConfig.containsKey('run_micro') ? (metadataPlotsConfig.run_micro as boolean) : true
    boolean metadataPlotsRunMito = metadataPlotsConfig.containsKey('run_mito') ? (metadataPlotsConfig.run_mito as boolean) : true
    if( metadataPlotsEnabled && !metadataPlotsRunMicro && !metadataPlotsRunMito ) {
        exit 1, "metadata_plots configured to skip both micro and mito outputs; disable metadata_plots.enabled instead."
    }
    boolean metadataForceMicroOnly = metadataPlotsRunMicro && !metadataPlotsRunMito
    boolean metadataForceMitoOnly = metadataPlotsRunMito && !metadataPlotsRunMicro

    def batchCorrectionConfig = config.batch_correction ?: [:]
    boolean batchCorrectionEnabled = metadataPlotsEnabled && (batchCorrectionConfig.containsKey('enabled') ? (batchCorrectionConfig.enabled as boolean) : true)
    def batchCorrectionOutputDir = batchCorrectionConfig.output_dir ?: 'batch_correction'
    def batchCorrectionOutputDirAbs = new File(outputDir, batchCorrectionOutputDir).canonicalPath
    def batchCorrectionBatchCol = batchCorrectionConfig.batch_col ?: 'batch'
    def batchCorrectionSampleIdCol = batchCorrectionConfig.sample_id_col ?: metadataPlotsSampleCol
    def batchCorrectionOrientation = batchCorrectionConfig.asv_orientation ?: 'features_rows'
    def batchBiologicalCovariates = batchCorrectionConfig.biological_covariates ? batchCorrectionConfig.biological_covariates.toString().trim() : ''
    def batchBiologicalColorCols = batchCorrectionConfig.biological_color_col ?: 'Depth'
    def batchColorPaletteCols = batchCorrectionConfig.color_palette_col ?: 'Color'
    def batchBiologicalPalettes = batchCorrectionConfig.biological_palettes instanceof Map ? batchCorrectionConfig.biological_palettes : [:]
    def batchBiologicalPalettesJson = groovy.json.JsonOutput.toJson(batchBiologicalPalettes)
    def batchUmapNeighbors = batchCorrectionConfig.umap_neighbors ? (batchCorrectionConfig.umap_neighbors as int) : 15
    def batchUmapMinDist = batchCorrectionConfig.umap_min_dist != null ? (batchCorrectionConfig.umap_min_dist as double) : 0.1d
    def batchHdbscanMinClusterSize = batchCorrectionConfig.hdbscan_min_cluster_size ? (batchCorrectionConfig.hdbscan_min_cluster_size as int) : 5
    def batchHdbscanMinSamples = batchCorrectionConfig.hdbscan_min_samples != null ? (batchCorrectionConfig.hdbscan_min_samples as int) : null
    def batchHdbscanSelectionMethod = batchCorrectionConfig.hdbscan_selection_method ?: 'eom'
    boolean batchOptimize = (batchCorrectionConfig.optimize_clustering ?: false) as boolean
    def batchTargetClusters = batchCorrectionConfig.target_clusters ?: '3-8'
    def batchNFeaturesPlot = batchCorrectionConfig.n_features_plot ? (batchCorrectionConfig.n_features_plot as int) : 5
    def batchRandomState = batchCorrectionConfig.random_state ? (batchCorrectionConfig.random_state as int) : 42
    def batchConqurMode = batchCorrectionConfig.conqur_mode ?: 'libsize'
    def batchConqurNumCore = batchCorrectionConfig.conqur_num_core ? (batchCorrectionConfig.conqur_num_core as int) : pipelineThreads
    def batchConqurBatchRef = batchCorrectionConfig.conqur_batch_ref ? batchCorrectionConfig.conqur_batch_ref.toString().trim() : ''
    boolean batchConqurLogisticLasso = (batchCorrectionConfig.conqur_logistic_lasso ?: false) as boolean
    def batchConqurQuantileType = batchCorrectionConfig.conqur_quantile_type ?: 'standard'
    boolean batchConqurSimpleMatch = (batchCorrectionConfig.conqur_simple_match ?: false) as boolean
    def batchConqurLambdaQuantile = batchCorrectionConfig.conqur_lambda_quantile ?: '2p/n'
    boolean batchConqurInterplt = (batchCorrectionConfig.conqur_interplt ?: false) as boolean
    def batchConqurDelta = batchCorrectionConfig.conqur_delta != null ? (batchCorrectionConfig.conqur_delta as double) : 0.4999d
    boolean batchConqurAutoInstall = (batchCorrectionConfig.conqur_auto_install ?: false) as boolean
    def batchCorrectionPolicy = batchCorrectionConfig.correction_policy ? batchCorrectionConfig.correction_policy.toString().trim().toLowerCase() : 'auto'
    if( !(batchCorrectionPolicy in ['auto','always','never']) ) {
        exit 1, "batch_correction.correction_policy must be one of: auto, always, never"
    }
    def batchAutoMinSampleRho = batchCorrectionConfig.auto_min_sample_rho != null ? (batchCorrectionConfig.auto_min_sample_rho as double) : 0.85d
    def batchAutoMinBrayRho = batchCorrectionConfig.auto_min_bray_rho != null ? (batchCorrectionConfig.auto_min_bray_rho as double) : 0.75d
    def batchAutoMaxBatchEtaRatio = batchCorrectionConfig.auto_max_batch_eta_ratio != null ? (batchCorrectionConfig.auto_max_batch_eta_ratio as double) : 0.95d
    def batchAutoMinBatchEtaDrop = batchCorrectionConfig.auto_min_batch_eta_drop != null ? (batchCorrectionConfig.auto_min_batch_eta_drop as double) : 0.01d
    def batchAutoMinBioEtaRatio = batchCorrectionConfig.auto_min_bio_eta_ratio != null ? (batchCorrectionConfig.auto_min_bio_eta_ratio as double) : 0.70d

    def outlierConfig = config.outlier_detection ?: [:]
    boolean outlierEnabled = batchCorrectionEnabled && (outlierConfig.containsKey('enabled') ? (outlierConfig.enabled as boolean) : true)
    def outlierOutputDir = outlierConfig.output_dir ?: 'outliers_corrected'
    def outlierOutputDirAbs = new File(outputDir, outlierOutputDir).canonicalPath
    def outlierSampleIdCol = outlierConfig.sample_id_col ?: (outlierConfig.sample_col ?: metadataPlotsSampleCol)
    def outlierGroupColsRaw = outlierConfig.group_cols
    List<String> outlierGroupCols = []
    if( outlierGroupColsRaw instanceof List ) {
        outlierGroupCols = outlierGroupColsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( outlierGroupColsRaw ) {
        outlierGroupCols = outlierGroupColsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    if( outlierGroupCols.isEmpty() ) {
        outlierGroupCols = ['none']
    }
    // OUTLIER_CHECKER consumes batch correction's CLR output (samples x ASVs).
    def outlierTransform = outlierConfig.transform ?: 'none'
    def outlierOrientation = outlierConfig.asv_orientation ?: 'samples_rows'
    boolean outlierPreTransformed = outlierConfig.containsKey('pre_transformed') ?
        (outlierConfig.pre_transformed as boolean) : true
    boolean outlierScale = (outlierConfig.scale ?: false) as boolean
    boolean outlierUseIso = (outlierConfig.use_iso ?: false) as boolean
    boolean outlierUseSvm = (outlierConfig.use_svm ?: false) as boolean
    boolean outlierUseHdb = (outlierConfig.use_hdb ?: false) as boolean
    def outlierVoteThreshold = outlierConfig.vote_threshold ? (outlierConfig.vote_threshold as int) : 3
    def outlierIsoContamination = outlierConfig.iso_contamination ?: 'auto'
    def outlierIsoEstimators = outlierConfig.iso_estimators ? (outlierConfig.iso_estimators as int) : 100
    def outlierIsoRandomState = outlierConfig.iso_random_state ? (outlierConfig.iso_random_state as int) : 42
    def outlierSvmKernel = outlierConfig.svm_kernel ?: 'rbf'
    def outlierSvmGamma = outlierConfig.svm_gamma ?: 'scale'
    def outlierSvmNu = outlierConfig.svm_nu ? (outlierConfig.svm_nu as double) : 0.1d
    def outlierHdbMinClusterSize = outlierConfig.hdbscan_min_cluster_size ? (outlierConfig.hdbscan_min_cluster_size as int) : 5
    def outlierHdbMinSamples = outlierConfig.hdbscan_min_samples != null ? (outlierConfig.hdbscan_min_samples as int) : null
    def outlierHdbMetric = outlierConfig.hdbscan_metric ?: 'euclidean'

    def collectorsConfig = config.collectors_curve ?: [:]
    boolean collectorsEnabled = metadataPlotsEnabled && (collectorsConfig.containsKey('enabled') ? (collectorsConfig.enabled as boolean) : true)
    def collectorsSampleCol = collectorsConfig.sample_col ?: metadataPlotsSampleCol
    def collectorsGroupCol = collectorsConfig.group_col ?: (collectorsConfig.group1_col ?: metadataPlotsTypeCol)
    def collectorsColorCol = collectorsConfig.color_col ?: 'Color'
    def collectorsGroupColors = collectorsConfig.group_colors ?: (collectorsConfig.group_palette ?: '')
    def collectorsGroupOrderRaw = collectorsConfig.group_order ?: metadataPlotsGroupOrder
    List<String> collectorsGroupOrder = normalizePresetList(collectorsGroupOrderRaw, [], config.order_presets ?: [:])
    def collectorsPermutations = collectorsConfig.permutations ? (collectorsConfig.permutations as int) : 999
    def collectorsSeed = collectorsConfig.seed ? (collectorsConfig.seed as int) : 42
    def collectorsOutPrefix = collectorsConfig.out_prefix ?: 'metadata/collectors_curve'
    def collectorsOutPrefixAbs = new File(outputDir, collectorsOutPrefix).canonicalPath
    def collectorsTitle = collectorsConfig.title ?: ''
    def collectorsFormats = collectorsConfig.formats ?: 'pdf,svg'
    def collectorsXpad = collectorsConfig.xpad != null ? (collectorsConfig.xpad as double) : 0.5d
    def collectorsMaxCols = collectorsConfig.max_cols ? (collectorsConfig.max_cols as int) : 3
    def collectorsShowPerms = collectorsConfig.show_perms ? (collectorsConfig.show_perms as int) : 10
    def collectorsPresenceThreshold = collectorsConfig.presence_threshold != null ? (collectorsConfig.presence_threshold as double) : 0d

    def plotUpsetConfig = config.plot_upset ?: [:]
    boolean plotUpsetRequested = plotUpsetConfig.containsKey('enabled') ? (plotUpsetConfig.enabled as boolean) : false
    if( plotUpsetRequested && !metadataPlotsEnabled ) {
        exit 1, "plot_upset.enabled requires metadata_plots.enabled to be true"
    }
    boolean plotUpsetEnabled = plotUpsetRequested
    def plotUpsetSubDir = plotUpsetConfig.sub_dir ?: '.'
    def plotUpsetDomain = plotUpsetConfig.domain ?: 'micro'
    def plotUpsetTaxonomyPath = plotUpsetConfig.taxonomy_path ? resolveOptionalPath(plotUpsetConfig.taxonomy_path, configRoot) : null
    def plotUpsetSampleIdCol = plotUpsetConfig.sample_id_col ?: metadataPlotsSampleCol
    def plotUpsetGroupCol = plotUpsetConfig.group_col ?: (plotUpsetConfig.group1_col ?: metadataPlotsTypeCol)
    def plotUpsetColorCol = plotUpsetConfig.color_col ?: metadataPlotsColorCol
    def plotUpsetGroupPalette = plotUpsetConfig.group_palette ?: ''
    def plotUpsetGroupOrderRaw = plotUpsetConfig.group_order ?: metadataPlotsGroupOrder
    List<String> plotUpsetGroupOrder = normalizePresetList(plotUpsetGroupOrderRaw, [], config.order_presets ?: [:])
    def plotUpsetSubsetGroupsRaw = plotUpsetConfig.subset_groups
    List<String> plotUpsetSubsetGroups = []
    if( plotUpsetSubsetGroupsRaw instanceof List ) {
        plotUpsetSubsetGroups = plotUpsetSubsetGroupsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( plotUpsetSubsetGroupsRaw ) {
        plotUpsetSubsetGroups = plotUpsetSubsetGroupsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def plotUpsetVennSubsetGroupsRaw = plotUpsetConfig.venn_subset_groups
    List<String> plotUpsetVennSubsetGroups = []
    if( plotUpsetVennSubsetGroupsRaw instanceof List ) {
        plotUpsetVennSubsetGroups = plotUpsetVennSubsetGroupsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( plotUpsetVennSubsetGroupsRaw ) {
        plotUpsetVennSubsetGroups = plotUpsetVennSubsetGroupsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    boolean plotUpsetSkipVenn = plotUpsetConfig.containsKey('skip_venn') ? (plotUpsetConfig.skip_venn as boolean) : true
    def plotUpsetMaxIntersections = plotUpsetConfig.max_intersections != null ? (plotUpsetConfig.max_intersections as int) : 20
    if( plotUpsetMaxIntersections < 1 ) {
        exit 1, "plot_upset.max_intersections must be at least 1"
    }
    def plotUpsetFormats = plotUpsetConfig.formats ?: 'pdf,svg,png'
    def plotUpsetFontSize = plotUpsetConfig.font_size != null ? (plotUpsetConfig.font_size as double) : 12d
    boolean plotUpsetRawOnly = plotUpsetConfig.containsKey('raw_only') ? (plotUpsetConfig.raw_only as boolean) : false
    boolean plotUpsetFinalOnly = plotUpsetConfig.containsKey('final_only') ? (plotUpsetConfig.final_only as boolean) : false
    if( plotUpsetRawOnly && plotUpsetFinalOnly ) {
        exit 1, "plot_upset.raw_only and plot_upset.final_only cannot both be true"
    }

    def bubbleplotterConfig = config.bubbleplotter ?: [:]
    boolean bubbleplotterRequested = bubbleplotterConfig.containsKey('enabled') ? (bubbleplotterConfig.enabled as boolean) : false
    if( bubbleplotterRequested && !metadataPlotsEnabled ) {
        exit 1, "bubbleplotter.enabled requires metadata_plots.enabled to be true"
    }
    boolean bubbleplotterEnabled = bubbleplotterRequested
    def bubbleplotterOutputPrefix = bubbleplotterConfig.output_prefix ?: 'metadata/bubble_plot_asv'
    def bubbleplotterOutputPrefixAbs = new File(outputDir, bubbleplotterOutputPrefix).canonicalPath
    def bubbleplotterOutputDirAbs = (new File(bubbleplotterOutputPrefixAbs).parentFile ?: new File(outputDir)).canonicalPath
    def bubbleplotterFormats = bubbleplotterConfig.formats ?: 'pdf,png,svg'
    def bubbleplotterCountCol = bubbleplotterConfig.count_col ?: 'count'
    def bubbleplotterSampleCol = bubbleplotterConfig.sample_col ?: metadataPlotsSampleCol
    def bubbleplotterDepthCol = bubbleplotterConfig.group1_col ?: (bubbleplotterConfig.depth_col ?: metadataPlotsTypeCol)
    def bubbleplotterColorCol = bubbleplotterConfig.color_col ?: metadataPlotsColorCol
    def bubbleplotterMonthCol = bubbleplotterConfig.group2_col ?: (bubbleplotterConfig.month_col ?: 'Month')
    def bubbleplotterGroup1OrderRaw = bubbleplotterConfig.group1_order ?: metadataPlotsGroupOrder
    List<String> bubbleplotterGroup1Order = normalizePresetList(bubbleplotterGroup1OrderRaw, [], config.order_presets ?: [:])
    def bubbleplotterGroup2OrderRaw = bubbleplotterConfig.group2_order ?: (config.indicspecies?.group2_order ?: '')
    List<String> bubbleplotterGroup2Order = normalizePresetList(bubbleplotterGroup2OrderRaw, [], config.order_presets ?: [:])
    def bubbleplotterFigsize = bubbleplotterConfig.figsize ?: '32,60'
    def bubbleplotterScale = bubbleplotterConfig.bubble_scale != null ? (bubbleplotterConfig.bubble_scale as double) : 10d
    boolean bubbleplotterNoAutoSize = bubbleplotterConfig.containsKey('no_auto_size') ? (bubbleplotterConfig.no_auto_size as boolean) : true

    def umapClusteringConfig = config.umap_clustering ?: [:]
    boolean umapClusteringRequested = umapClusteringConfig.containsKey('enabled') ? (umapClusteringConfig.enabled as boolean) : false
    if( umapClusteringRequested && !metadataPlotsEnabled ) {
        exit 1, "umap_clustering.enabled requires metadata_plots.enabled to be true"
    }
    boolean umapClusteringEnabled = umapClusteringRequested
    def umapClusteringOutputPrefix = umapClusteringConfig.output_prefix ?: 'metadata/umap_clustering'
    def umapClusteringOutputPrefixAbs = new File(outputDir, umapClusteringOutputPrefix).canonicalPath
    def umapClusteringOutputDirAbs = (new File(umapClusteringOutputPrefixAbs).parentFile ?: new File(outputDir)).canonicalPath
    def umapClusteringSampleCol = umapClusteringConfig.sample_col ?: metadataPlotsSampleCol
    def umapClusteringCountCol = umapClusteringConfig.count_col ?: 'count'
    def umapClusteringDepthCol = umapClusteringConfig.group1_col ?: (umapClusteringConfig.depth_col ?: metadataPlotsTypeCol)
    def umapClusteringColorCol = umapClusteringConfig.color_col ?: metadataPlotsColorCol
    def umapClusteringSecondaryCol = umapClusteringConfig.group2_col ?: (umapClusteringConfig.secondary_col ?: (umapClusteringConfig.month_col ?: 'Month'))
    def umapClusteringGroup1OrderRaw = umapClusteringConfig.group1_order ?: metadataPlotsGroupOrder
    List<String> umapClusteringGroup1Order = normalizePresetList(umapClusteringGroup1OrderRaw, [], config.order_presets ?: [:])
    def umapClusteringGroup2OrderRaw = umapClusteringConfig.group2_order ?: (config.indicspecies?.group2_order ?: '')
    List<String> umapClusteringGroup2Order = normalizePresetList(umapClusteringGroup2OrderRaw, [], config.order_presets ?: [:])
    def sharedPaletteConfig = config.indicspecies?.group_palettes instanceof Map ? config.indicspecies.group_palettes : [:]
    def umapClusteringGroup1Palette = umapClusteringConfig.group1_palette ?: (sharedPaletteConfig[umapClusteringDepthCol] ?: '')
    def umapClusteringGroup2Palette = umapClusteringConfig.group2_palette ?: (sharedPaletteConfig[umapClusteringSecondaryCol] ?: '')
    def umapClusteringFormats = umapClusteringConfig.formats ?: 'pdf,png,svg'
    def umapClusteringNormalize = umapClusteringConfig.normalize ?: 'clr'
    def umapClusteringTransform = umapClusteringConfig.transform ?: 'sqrt'
    def umapClusteringNeighbors = umapClusteringConfig.n_neighbors ? (umapClusteringConfig.n_neighbors as int) : 15
    def umapClusteringMinDist = umapClusteringConfig.min_dist != null ? (umapClusteringConfig.min_dist as double) : 0.1d
    def umapClusteringMetric = umapClusteringConfig.umap_metric ?: 'euclidean'
    def umapClusteringHdbscanMetric = umapClusteringConfig.hdbscan_metric ?: 'euclidean'
    def umapClusteringMinClusterSize = umapClusteringConfig.min_cluster_size ? (umapClusteringConfig.min_cluster_size as int) : 10
    def umapClusteringMinSamples = umapClusteringConfig.min_samples ? (umapClusteringConfig.min_samples as int) : 5
    boolean umapClusteringNoScale = (umapClusteringConfig.no_scale ?: false) as boolean
    def umapClusteringRandomState = umapClusteringConfig.random_state ? (umapClusteringConfig.random_state as int) : 42

    def diversityConfig = config.diversity ?: [:]
    boolean diversityRequested = diversityConfig.containsKey('enabled') ? (diversityConfig.enabled as boolean) : false
    if( diversityRequested && !metadataPlotsEnabled ) {
        exit 1, "diversity.enabled requires metadata_plots.enabled to be true"
    }
    boolean diversityEnabled = diversityRequested
    def diversityOutputDir = diversityConfig.output_dir ?: 'diversity'
    def diversityOutputDirAbs = new File(outputDir, diversityOutputDir).canonicalPath
    def diversityMatchedConfig = diversityConfig.matched_cohort ?: [:]
    def diversityMatchedCohortPath = diversityMatchedConfig.table ? resolveOptionalPath(diversityMatchedConfig.table, configRoot) : ''
    def diversityMatchedCohortSource = diversityMatchedConfig.source ?: 'external'
    if( !(diversityMatchedCohortSource in ['external', 'community_predictor_comparison']) ) {
        exit 1, "diversity.matched_cohort.source must be external or community_predictor_comparison"
    }
    if( diversityMatchedCohortSource == 'community_predictor_comparison' && diversityMatchedCohortPath ) {
        exit 1, "Set either a generated diversity matched-cohort source or an external table, not both"
    }
    def diversityMitoOutputDir = diversityConfig.mito_output_dir ?: 'mito/diversity'
    def diversityMitoOutputDirAbs = new File(outputDir, diversityMitoOutputDir).canonicalPath
    def diversityMitoInputPath = diversityConfig.mito_input ? resolveOptionalPath(diversityConfig.mito_input, configRoot) : new File(outputDir, 'mito/ASVs/ASV_target.mito.tsv').canonicalPath
    def diversitySampleCol = diversityConfig.sample_col ?: metadataPlotsSampleCol
    def diversityGroupCol = diversityConfig.group_col ?: (diversityConfig.group1_col ?: metadataPlotsTypeCol)
    def diversityColorCol = diversityConfig.color_col ?: 'Color'
    def diversitySecondaryCol = diversityConfig.group2_col ?: (diversityConfig.secondary_col ?: 'Month')
    def diversityGroupPalette = diversityConfig.group1_palette ?: (sharedPaletteConfig[diversityGroupCol] ?: '')
    def diversitySecondaryPalette = diversityConfig.group2_palette ?: (sharedPaletteConfig[diversitySecondaryCol] ?: '')
    def diversityExcludeGroupsRaw = diversityConfig.exclude_groups
    List<String> diversityExcludeGroups = []
    if( diversityExcludeGroupsRaw instanceof List ) {
        diversityExcludeGroups = diversityExcludeGroupsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( diversityExcludeGroupsRaw ) {
        diversityExcludeGroups = diversityExcludeGroupsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def diversityGroupOrderRaw = diversityConfig.group_order ?: metadataPlotsGroupOrder
    List<String> diversityGroupOrder = normalizePresetList(diversityGroupOrderRaw, [], config.order_presets ?: [:])
    boolean diversityRunMito = diversityConfig.containsKey('run_mito') ? (diversityConfig.run_mito as boolean) : true
    def diversityUmapNeighbors = diversityConfig.umap_neighbors ? (diversityConfig.umap_neighbors as int) : 30
    def diversityUmapMinDist = diversityConfig.umap_min_dist != null ? (diversityConfig.umap_min_dist as double) : 0.01d
    def diversityPermutations = diversityConfig.permanova_perms ? (diversityConfig.permanova_perms as int) : 999
    def diversityRandomState = diversityConfig.random_state ? (diversityConfig.random_state as int) : 42
    def diversityBlockCol = diversityConfig.block_col ? diversityConfig.block_col.toString().trim() : ''
    boolean diversityVerbose = diversityConfig.containsKey('verbose') ? (diversityConfig.verbose as boolean) : true
    def diversityPatientAwareConfig = (diversityConfig.patient_aware instanceof Map) ? diversityConfig.patient_aware : [:]
    boolean diversityPatientAwareEnabled = diversityPatientAwareConfig.containsKey('enabled') ? (diversityPatientAwareConfig.enabled as boolean) : false
    def diversityPatientAwareOutputDir = diversityPatientAwareConfig.output_dir ?: 'patient_aware'
    def diversityPatientAwareOutputDirAbs = new File(diversityOutputDirAbs, diversityPatientAwareOutputDir.toString()).canonicalPath
    def diversityPatientAwareSampleCol = diversityPatientAwareConfig.sample_col ?: diversitySampleCol
    def diversityPatientAwarePatientCol = diversityPatientAwareConfig.patient_col ?
        diversityPatientAwareConfig.patient_col.toString().trim() :
        (diversityBlockCol ?: 'Participant_ID')
    def diversityPatientAwareCaseCol = diversityPatientAwareConfig.case_col ?
        diversityPatientAwareConfig.case_col.toString().trim() : 'Case'
    def diversityPatientAwareTypeCol = diversityPatientAwareConfig.type_col ?: diversityGroupCol
    def diversityPatientAwareSampleTypesRaw = diversityPatientAwareConfig.sample_types ?:
        (diversityGroupOrder ? diversityGroupOrder.join(',') : 'Oral Rinse,BAL,Lung Brush')
    def diversityPatientAwareSampleTypes = diversityPatientAwareSampleTypesRaw instanceof List ?
        diversityPatientAwareSampleTypesRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        diversityPatientAwareSampleTypesRaw.toString().trim()
    boolean diversityPatientAwareExcludeContralateral = diversityPatientAwareConfig.containsKey('exclude_contralateral_in_cancer') ?
        (diversityPatientAwareConfig.exclude_contralateral_in_cancer as boolean) : true
    def diversityPatientAwareContralateralCol = diversityPatientAwareConfig.contralateral_col ?
        diversityPatientAwareConfig.contralateral_col.toString().trim() : 'lung_status'
    def diversityPatientAwareCancerSiteCol = diversityPatientAwareConfig.cancer_site_col ?
        diversityPatientAwareConfig.cancer_site_col.toString().trim() : 'Cancer_Site'
    def diversityPatientAwareLungSideCol = diversityPatientAwareConfig.lung_side_col ?
        diversityPatientAwareConfig.lung_side_col.toString().trim() : 'lung_code'
    def diversityPatientAwareContralateralValue = diversityPatientAwareConfig.contralateral_value ?
        diversityPatientAwareConfig.contralateral_value.toString().trim() : 'Contralateral'
    def diversityPatientAwareContralateralTypesRaw = diversityPatientAwareConfig.contralateral_sample_types ?: 'Lung Brush,BAL'
    def diversityPatientAwareContralateralTypes = diversityPatientAwareContralateralTypesRaw instanceof List ?
        diversityPatientAwareContralateralTypesRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        diversityPatientAwareContralateralTypesRaw.toString().trim()
    def diversityPatientAwareTransform = diversityPatientAwareConfig.transform ?
        diversityPatientAwareConfig.transform.toString().trim().toLowerCase() : 'none'
    if( !['none', 'rclr'].contains(diversityPatientAwareTransform) ) {
        exit 1, "diversity.patient_aware.transform must be one of: none, rclr"
    }
    def diversityPatientAwarePermutations = diversityPatientAwareConfig.permutations ?
        (diversityPatientAwareConfig.permutations as int) : 9999
    def diversityPatientAwareSeed = diversityPatientAwareConfig.seed ?
        (diversityPatientAwareConfig.seed as int) : diversityRandomState
    boolean diversityPatientAwareRequireCompleteTypes = diversityPatientAwareConfig.containsKey('require_complete_types') ?
        (diversityPatientAwareConfig.require_complete_types as boolean) : false
    return [
        metadataPlotsMetadataPath: metadataPlotsMetadataPath,
        metadataPlotsSubDir: metadataPlotsSubDir,
        metadataPlotsSampleCol: metadataPlotsSampleCol,
        metadataPlotsTypeCol: metadataPlotsTypeCol,
        metadataPlotsColorCol: metadataPlotsColorCol,
        metadataPlotsBiochemAssignmentsPath: metadataPlotsBiochemAssignmentsPath,
        metadataPlotsBiochemSampleCol: metadataPlotsBiochemSampleCol,
        metadataPlotsStratificationTimeseriesPath: metadataPlotsStratificationTimeseriesPath,
        metadataPlotsStratMetaJoinCol: metadataPlotsStratMetaJoinCol,
        metadataPlotsStratJoinCol: metadataPlotsStratJoinCol,
        metadataPlotsCruiseGroupPath: metadataPlotsCruiseGroupPath,
        metadataPlotsCruiseGroupMetaJoinCol: metadataPlotsCruiseGroupMetaJoinCol,
        metadataPlotsCruiseGroupJoinCol: metadataPlotsCruiseGroupJoinCol,
        metadataPlotsSubtractionCol: metadataPlotsSubtractionCol,
        batchCorrectionOutputDir: batchCorrectionOutputDir,
        batchCorrectionOutputDirAbs: batchCorrectionOutputDirAbs,
        batchCorrectionBatchCol: batchCorrectionBatchCol,
        batchCorrectionSampleIdCol: batchCorrectionSampleIdCol,
        batchCorrectionOrientation: batchCorrectionOrientation,
        batchBiologicalCovariates: batchBiologicalCovariates,
        batchBiologicalColorCols: batchBiologicalColorCols,
        batchColorPaletteCols: batchColorPaletteCols,
        batchBiologicalPalettesJson: batchBiologicalPalettesJson,
        batchUmapNeighbors: batchUmapNeighbors,
        batchUmapMinDist: batchUmapMinDist,
        batchHdbscanMinClusterSize: batchHdbscanMinClusterSize,
        batchHdbscanMinSamples: batchHdbscanMinSamples,
        batchHdbscanSelectionMethod: batchHdbscanSelectionMethod,
        batchTargetClusters: batchTargetClusters,
        batchNFeaturesPlot: batchNFeaturesPlot,
        batchRandomState: batchRandomState,
        batchConqurMode: batchConqurMode,
        batchConqurNumCore: batchConqurNumCore,
        batchConqurBatchRef: batchConqurBatchRef,
        batchConqurQuantileType: batchConqurQuantileType,
        batchConqurLambdaQuantile: batchConqurLambdaQuantile,
        batchConqurDelta: batchConqurDelta,
        batchCorrectionPolicy: batchCorrectionPolicy,
        batchAutoMinSampleRho: batchAutoMinSampleRho,
        batchAutoMinBrayRho: batchAutoMinBrayRho,
        batchAutoMaxBatchEtaRatio: batchAutoMaxBatchEtaRatio,
        batchAutoMinBatchEtaDrop: batchAutoMinBatchEtaDrop,
        batchAutoMinBioEtaRatio: batchAutoMinBioEtaRatio,
        outlierOutputDirAbs: outlierOutputDirAbs,
        outlierSampleIdCol: outlierSampleIdCol,
        outlierTransform: outlierTransform,
        outlierOrientation: outlierOrientation,
        outlierVoteThreshold: outlierVoteThreshold,
        outlierIsoContamination: outlierIsoContamination,
        outlierIsoEstimators: outlierIsoEstimators,
        outlierIsoRandomState: outlierIsoRandomState,
        outlierSvmKernel: outlierSvmKernel,
        outlierSvmGamma: outlierSvmGamma,
        outlierSvmNu: outlierSvmNu,
        outlierHdbMinClusterSize: outlierHdbMinClusterSize,
        outlierHdbMinSamples: outlierHdbMinSamples,
        outlierHdbMetric: outlierHdbMetric,
        collectorsSampleCol: collectorsSampleCol,
        collectorsGroupCol: collectorsGroupCol,
        collectorsColorCol: collectorsColorCol,
        collectorsPermutations: collectorsPermutations,
        collectorsSeed: collectorsSeed,
        collectorsOutPrefixAbs: collectorsOutPrefixAbs,
        collectorsTitle: collectorsTitle,
        collectorsFormats: collectorsFormats,
        collectorsXpad: collectorsXpad,
        collectorsMaxCols: collectorsMaxCols,
        collectorsShowPerms: collectorsShowPerms,
        collectorsPresenceThreshold: collectorsPresenceThreshold,
        plotUpsetSubDir: plotUpsetSubDir,
        plotUpsetDomain: plotUpsetDomain,
        plotUpsetTaxonomyPath: plotUpsetTaxonomyPath,
        plotUpsetSampleIdCol: plotUpsetSampleIdCol,
        plotUpsetGroupCol: plotUpsetGroupCol,
        plotUpsetColorCol: plotUpsetColorCol,
        plotUpsetFormats: plotUpsetFormats,
        plotUpsetFontSize: plotUpsetFontSize,
        bubbleplotterOutputPrefixAbs: bubbleplotterOutputPrefixAbs,
        bubbleplotterOutputDirAbs: bubbleplotterOutputDirAbs,
        bubbleplotterFormats: bubbleplotterFormats,
        bubbleplotterCountCol: bubbleplotterCountCol,
        bubbleplotterSampleCol: bubbleplotterSampleCol,
        bubbleplotterDepthCol: bubbleplotterDepthCol,
        bubbleplotterColorCol: bubbleplotterColorCol,
        bubbleplotterMonthCol: bubbleplotterMonthCol,
        bubbleplotterFigsize: bubbleplotterFigsize,
        bubbleplotterScale: bubbleplotterScale,
        umapClusteringOutputPrefixAbs: umapClusteringOutputPrefixAbs,
        umapClusteringOutputDirAbs: umapClusteringOutputDirAbs,
        umapClusteringSampleCol: umapClusteringSampleCol,
        umapClusteringCountCol: umapClusteringCountCol,
        umapClusteringDepthCol: umapClusteringDepthCol,
        umapClusteringColorCol: umapClusteringColorCol,
        umapClusteringSecondaryCol: umapClusteringSecondaryCol,
        umapClusteringGroup1Palette: umapClusteringGroup1Palette,
        umapClusteringGroup2Palette: umapClusteringGroup2Palette,
        umapClusteringFormats: umapClusteringFormats,
        umapClusteringNormalize: umapClusteringNormalize,
        umapClusteringTransform: umapClusteringTransform,
        umapClusteringNeighbors: umapClusteringNeighbors,
        umapClusteringMinDist: umapClusteringMinDist,
        umapClusteringMetric: umapClusteringMetric,
        umapClusteringHdbscanMetric: umapClusteringHdbscanMetric,
        umapClusteringMinClusterSize: umapClusteringMinClusterSize,
        umapClusteringMinSamples: umapClusteringMinSamples,
        umapClusteringRandomState: umapClusteringRandomState,
        diversityOutputDirAbs: diversityOutputDirAbs,
        diversityMatchedCohortPath: diversityMatchedCohortPath,
        diversityMatchedCohortSource: diversityMatchedCohortSource,
        diversityMitoOutputDirAbs: diversityMitoOutputDirAbs,
        diversityMitoInputPath: diversityMitoInputPath,
        diversitySampleCol: diversitySampleCol,
        diversityGroupCol: diversityGroupCol,
        diversityColorCol: diversityColorCol,
        diversitySecondaryCol: diversitySecondaryCol,
        diversityGroupPalette: diversityGroupPalette,
        diversitySecondaryPalette: diversitySecondaryPalette,
        diversityUmapNeighbors: diversityUmapNeighbors,
        diversityUmapMinDist: diversityUmapMinDist,
        diversityPermutations: diversityPermutations,
        diversityRandomState: diversityRandomState,
        diversityBlockCol: diversityBlockCol,
        diversityPatientAwareOutputDirAbs: diversityPatientAwareOutputDirAbs,
        diversityPatientAwareSampleCol: diversityPatientAwareSampleCol,
        diversityPatientAwarePatientCol: diversityPatientAwarePatientCol,
        diversityPatientAwareCaseCol: diversityPatientAwareCaseCol,
        diversityPatientAwareTypeCol: diversityPatientAwareTypeCol,
        diversityPatientAwareSampleTypes: diversityPatientAwareSampleTypes,
        diversityPatientAwareContralateralCol: diversityPatientAwareContralateralCol,
        diversityPatientAwareCancerSiteCol: diversityPatientAwareCancerSiteCol,
        diversityPatientAwareLungSideCol: diversityPatientAwareLungSideCol,
        diversityPatientAwareContralateralValue: diversityPatientAwareContralateralValue,
        diversityPatientAwareContralateralTypes: diversityPatientAwareContralateralTypes,
        diversityPatientAwareTransform: diversityPatientAwareTransform,
        diversityPatientAwarePermutations: diversityPatientAwarePermutations,
        diversityPatientAwareSeed: diversityPatientAwareSeed,
        metadataPlotsEnabled: metadataPlotsEnabled,
        metadataPlotsBiochemIncludeCols: metadataPlotsBiochemIncludeCols,
        metadataPlotsStratIncludeCols: metadataPlotsStratIncludeCols,
        metadataPlotsCruiseGroupIncludeCols: metadataPlotsCruiseGroupIncludeCols,
        metadataPlotsBiochemMetaJoinCols: metadataPlotsBiochemMetaJoinCols,
        metadataPlotsBiochemJoinCols: metadataPlotsBiochemJoinCols,
        metadataGroupNormalizationEnabled: metadataGroupNormalizationEnabled,
        metadataGroupNormalizationCols: metadataGroupNormalizationCols,
        metadataGroupNormalizationPattern: metadataGroupNormalizationPattern,
        metadataGroupNormalizationReplacement: metadataGroupNormalizationReplacement,
        metadataGroupNormalizationPreserveSource: metadataGroupNormalizationPreserveSource,
        metadataKeepTypes: metadataKeepTypes,
        metadataPlotsSubtractionGroups: metadataPlotsSubtractionGroups,
        metadataPlotsGroupOrder: metadataPlotsGroupOrder,
        metadataIncludeRank: metadataIncludeRank,
        batchCorrectionEnabled: batchCorrectionEnabled,
        batchOptimize: batchOptimize,
        batchConqurLogisticLasso: batchConqurLogisticLasso,
        batchConqurSimpleMatch: batchConqurSimpleMatch,
        batchConqurInterplt: batchConqurInterplt,
        batchConqurAutoInstall: batchConqurAutoInstall,
        outlierEnabled: outlierEnabled,
        outlierGroupCols: outlierGroupCols,
        outlierPreTransformed: outlierPreTransformed,
        outlierScale: outlierScale,
        outlierUseIso: outlierUseIso,
        outlierUseSvm: outlierUseSvm,
        outlierUseHdb: outlierUseHdb,
        collectorsEnabled: collectorsEnabled,
        collectorsGroupColors: collectorsGroupColors,
        collectorsGroupOrder: collectorsGroupOrder,
        plotUpsetEnabled: plotUpsetEnabled,
        plotUpsetGroupPalette: plotUpsetGroupPalette,
        plotUpsetGroupOrder: plotUpsetGroupOrder,
        plotUpsetSubsetGroups: plotUpsetSubsetGroups,
        plotUpsetVennSubsetGroups: plotUpsetVennSubsetGroups,
        plotUpsetSkipVenn: plotUpsetSkipVenn,
        plotUpsetMaxIntersections: plotUpsetMaxIntersections,
        plotUpsetRawOnly: plotUpsetRawOnly,
        plotUpsetFinalOnly: plotUpsetFinalOnly,
        bubbleplotterEnabled: bubbleplotterEnabled,
        bubbleplotterGroup1Order: bubbleplotterGroup1Order,
        bubbleplotterGroup2Order: bubbleplotterGroup2Order,
        bubbleplotterNoAutoSize: bubbleplotterNoAutoSize,
        umapClusteringEnabled: umapClusteringEnabled,
        umapClusteringGroup1Order: umapClusteringGroup1Order,
        umapClusteringGroup2Order: umapClusteringGroup2Order,
        umapClusteringNoScale: umapClusteringNoScale,
        diversityEnabled: diversityEnabled,
        diversityExcludeGroups: diversityExcludeGroups,
        diversityGroupOrder: diversityGroupOrder,
        diversityRunMito: diversityRunMito,
        diversityVerbose: diversityVerbose,
        diversityPatientAwareEnabled: diversityPatientAwareEnabled,
        diversityPatientAwareExcludeContralateral: diversityPatientAwareExcludeContralateral,
        diversityPatientAwareRequireCompleteTypes: diversityPatientAwareRequireCompleteTypes,
    ]
}

def parseIndicatorAndNetworkConfig(config, File configRoot, String outputDir, int pipelineThreads, Map analysisCfg) {
    def metadataPlotsEnabled = analysisCfg.metadataPlotsEnabled
    def metadataPlotsMetadataPath = analysisCfg.metadataPlotsMetadataPath
    def metadataPlotsSampleCol = analysisCfg.metadataPlotsSampleCol
    def metadataPlotsTypeCol = analysisCfg.metadataPlotsTypeCol
    def metadataPlotsColorCol = analysisCfg.metadataPlotsColorCol
    def metadataPlotsGroupOrder = analysisCfg.metadataPlotsGroupOrder
    def indicspeciesConfig = config.indicspecies ?: [:]
    boolean indicspeciesRequested = indicspeciesConfig.containsKey('enabled') ? (indicspeciesConfig.enabled as boolean) : false
    if( indicspeciesRequested && !metadataPlotsEnabled ) {
        exit 1, "indicspecies.enabled requires metadata_plots.enabled to be true"
    }
    def indicspeciesGroupColsRaw = indicspeciesConfig.group_cols
    List<String> indicspeciesGroupCols = []
    if( indicspeciesGroupColsRaw instanceof List ) {
        indicspeciesGroupCols = indicspeciesGroupColsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( indicspeciesGroupColsRaw ) {
        indicspeciesGroupCols = indicspeciesGroupColsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    } else {
        indicspeciesGroupCols = ['Depth', 'Month']
    }
    if( indicspeciesRequested && indicspeciesGroupCols.size() < 2 ) {
        exit 1, "indicspecies.group_cols must contain at least two groups when indicspecies.enabled is true"
    }
    boolean indicspeciesEnabled = indicspeciesRequested
    def indicspeciesSampleCol = indicspeciesConfig.sample_col ?: metadataPlotsSampleCol
    def indicspeciesPerms = indicspeciesConfig.perms ? (indicspeciesConfig.perms as int) : 9999
    def indicspeciesSeed = indicspeciesConfig.seed ? (indicspeciesConfig.seed as int) : 42
    def indicspeciesQThreshold = indicspeciesConfig.q_threshold != null ? (indicspeciesConfig.q_threshold as double) : 0.05d
    if( indicspeciesQThreshold < 0 || indicspeciesQThreshold > 1 ) {
        exit 1, "indicspecies.q_threshold must be between 0 and 1"
    }
    def indicspeciesMinN = indicspeciesConfig.min_n ? (indicspeciesConfig.min_n as int) : 2
    def indicspeciesBlockCol = indicspeciesConfig.block_col ? indicspeciesConfig.block_col.toString().trim() : ''
    def indicspeciesBlockedColsRaw = indicspeciesConfig.blocked_cols ?: []
    def indicspeciesBlockedCols = indicspeciesBlockedColsRaw instanceof List ?
        indicspeciesBlockedColsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        indicspeciesBlockedColsRaw.toString().trim()
    def indicspeciesStratifiedConfig = indicspeciesConfig.stratified ?: [:]
    boolean indicspeciesStratifiedEnabled = indicspeciesStratifiedConfig instanceof Map ?
        (indicspeciesStratifiedConfig.containsKey('enabled') ? (indicspeciesStratifiedConfig.enabled as boolean) : false) :
        false
    List<String> indicspeciesStratifiedSpecs = []
    if( indicspeciesStratifiedEnabled ) {
        def stratifiedAnalysesRaw = indicspeciesStratifiedConfig.analyses ?: []
        if( !(stratifiedAnalysesRaw instanceof List) ) {
            exit 1, "indicspecies.stratified.analyses must be a list when indicspecies.stratified.enabled is true"
        }
        indicspeciesStratifiedSpecs = stratifiedAnalysesRaw.collect { analysis ->
            if( !(analysis instanceof Map) ) {
                exit 1, "Each indicspecies.stratified.analyses entry must be a map with within_col and group_col"
            }
            def withinCol = (analysis.within_col ?: analysis.within ?: '').toString().trim()
            def groupCol = (analysis.group_col ?: analysis.group ?: '').toString().trim()
            if( !withinCol || !groupCol ) {
                exit 1, "Each indicspecies.stratified.analyses entry requires within_col and group_col"
            }
            def unsafeValues = [withinCol, groupCol].findAll { it.contains('::') || it.contains(';') || it.contains('|') }
            if( unsafeValues ) {
                exit 1, "indicspecies.stratified column names cannot contain '::', ';', or '|': ${unsafeValues.join(', ')}"
            }
            def levelsRaw = analysis.levels ?: analysis.within_values ?: analysis.sample_types ?: []
            List<String> levels = []
            if( levelsRaw instanceof List ) {
                levels = levelsRaw.collect { it.toString().trim() }.findAll { it }
            } else if( levelsRaw ) {
                levels = levelsRaw.toString().split(/\|/).collect { it.trim() }.findAll { it }
            }
            def unsafeLevels = levels.findAll { it.contains('::') || it.contains(';') || it.contains('|') }
            if( unsafeLevels ) {
                exit 1, "indicspecies.stratified levels cannot contain '::', ';', or '|': ${unsafeLevels.join(', ')}"
            }
            levels ? "${withinCol}::${groupCol}::${levels.join('|')}" : "${withinCol}::${groupCol}"
        }.findAll { it }
        if( indicspeciesStratifiedSpecs.isEmpty() ) {
            exit 1, "indicspecies.stratified.enabled is true but no valid analyses were configured"
        }
    }
    def indicspeciesStratifiedSpecsArg = indicspeciesStratifiedSpecs.join(';')
    def indicspeciesGroup1 = indicspeciesGroupCols ? indicspeciesGroupCols[0] : 'group1'
    def indicspeciesGroup2 = indicspeciesGroupCols.size() > 1 ? indicspeciesGroupCols[1] : ''
    def indicspeciesOutputDirAbs = new File(outputDir, "indicspecies").canonicalPath
    boolean indicspeciesPlotEnabled = indicspeciesConfig.containsKey('plot_enabled') ? (indicspeciesConfig.plot_enabled as boolean) : true
    def indicspeciesPlotPairsMode = indicspeciesConfig.plot_pairs_mode ?: 'all'
    def indicspeciesPlotOutputDir = indicspeciesConfig.plot_output_dir ?: 'indicspecies/plots'
    def indicspeciesPlotOutputDirAbs = new File(outputDir, indicspeciesPlotOutputDir).canonicalPath
    def indicspeciesPlotVennPath = indicspeciesConfig.venn ? resolveOptionalPath(indicspeciesConfig.venn, configRoot) : null
    def indicspeciesPlotTaxonomyPath = indicspeciesConfig.taxonomy ? resolveOptionalPath(indicspeciesConfig.taxonomy, configRoot) : new File(outputDir, 'taxonomy/ASV_SILVA_tax.full-length.vsearch.tsv').canonicalPath
    def indicspeciesColorCol = indicspeciesConfig.color_col ?: metadataPlotsColorCol
    def indicspeciesGroupPaletteMap = extractNamedStringMap(indicspeciesConfig as Map, indicspeciesGroupCols, 'group_palettes', 'palette')
    def indicspeciesGroupOrderMap = extractNamedListMap(indicspeciesConfig as Map, indicspeciesGroupCols, 'group_orders', 'order')
    def indicspeciesFocusLabelMap = extractNamedStringMap(indicspeciesConfig as Map, indicspeciesGroupCols, 'focus_labels', 'focus_label')
    if( indicspeciesGroup1 && !indicspeciesGroupOrderMap.containsKey(indicspeciesGroup1) && metadataPlotsGroupOrder ) {
        indicspeciesGroupOrderMap[indicspeciesGroup1] = metadataPlotsGroupOrder
    }
    def indicspeciesGroup1Palette = indicspeciesGroupPaletteMap[indicspeciesGroup1] ?: ''
    def indicspeciesGroup2Palette = indicspeciesGroup2 ? (indicspeciesGroupPaletteMap[indicspeciesGroup2] ?: '') : ''
    List<String> indicspeciesGroup1Order = indicspeciesGroupOrderMap[indicspeciesGroup1] ?: []
    List<String> indicspeciesGroup2Order = indicspeciesGroup2 ? (indicspeciesGroupOrderMap[indicspeciesGroup2] ?: []) : []
    def indicspeciesFocusGroup1Label = indicspeciesFocusLabelMap[indicspeciesGroup1] ?: ''
    def indicspeciesFocusGroup2Label = indicspeciesGroup2 ? (indicspeciesFocusLabelMap[indicspeciesGroup2] ?: '') : ''
    def indicspeciesGroupColsCsv = indicspeciesGroupCols.join(',')
    def indicspeciesGroupPaletteJson = groovy.json.JsonOutput.toJson(indicspeciesGroupPaletteMap)
    def indicspeciesGroupOrderJson = groovy.json.JsonOutput.toJson(indicspeciesGroupOrderMap)
    def indicspeciesFocusLabelJson = groovy.json.JsonOutput.toJson(indicspeciesFocusLabelMap)
    boolean indicspeciesLabelFocusedAsvs = indicspeciesConfig.containsKey('label_focused_asvs') ? (indicspeciesConfig.label_focused_asvs as boolean) : false
    boolean indicspeciesAlignedEnabled = indicspeciesConfig.containsKey('aligned_plot_enabled') ? (indicspeciesConfig.aligned_plot_enabled as boolean) : false
    def indicspeciesAlignedOutputDir = indicspeciesConfig.aligned_plot_output_dir ?: 'indicspecies/aligned'
    def indicspeciesAlignedOutputDirAbs = new File(outputDir, indicspeciesAlignedOutputDir).canonicalPath
    def indicspeciesAlignedAlpha = indicspeciesConfig.aligned_alpha != null ? (indicspeciesConfig.aligned_alpha as double) : 0.05d
    def indicspeciesAlignedMinStat = indicspeciesConfig.aligned_min_stat != null ? (indicspeciesConfig.aligned_min_stat as double) : 0.0d
    def indicspeciesAlignedTopN = indicspeciesConfig.aligned_top_n ? (indicspeciesConfig.aligned_top_n as int) : 25

    def vocCorrelationConfig = config.voc_correlation ?: [:]
    boolean vocCorrelationRequested = vocCorrelationConfig.containsKey('enabled') ? (vocCorrelationConfig.enabled as boolean) : false
    if( vocCorrelationRequested && !metadataPlotsEnabled ) {
        exit 1, "voc_correlation.enabled requires metadata_plots.enabled to be true"
    }
    def vocCorrelationVocTablePath = vocCorrelationConfig.voc_table ? resolveOptionalPath(vocCorrelationConfig.voc_table, configRoot) : null
    if( vocCorrelationRequested && (!vocCorrelationVocTablePath || !new File(vocCorrelationVocTablePath).exists()) ) {
        exit 1, "voc_correlation.voc_table must point to an existing VOC table when voc_correlation.enabled is true"
    }
    boolean vocCorrelationEnabled = vocCorrelationRequested
    def vocCorrelationOutputDir = vocCorrelationConfig.output_dir ?: 'voc_correlation'
    def vocCorrelationOutputDirAbs = resolveOutputRelative(vocCorrelationOutputDir.toString(), outputDir)
    def vocCorrelationVocSampleCol = vocCorrelationConfig.voc_sample_col ? vocCorrelationConfig.voc_sample_col.toString().trim() : 'sample'
    def vocCorrelationSampleIdMode = vocCorrelationConfig.sample_id_mode ? vocCorrelationConfig.sample_id_mode.toString().trim() : 'legacy_patient_pair'
    def vocCorrelationMetadataSampleCol = vocCorrelationConfig.metadata_sample_col ? vocCorrelationConfig.metadata_sample_col.toString().trim() : metadataPlotsSampleCol
    def vocCorrelationTypeCol = vocCorrelationConfig.type_col ? vocCorrelationConfig.type_col.toString().trim() : metadataPlotsTypeCol
    def vocCorrelationPatientCol = vocCorrelationConfig.patient_col ? vocCorrelationConfig.patient_col.toString().trim() : 'Participant_ID'
    def vocCorrelationCaseCol = vocCorrelationConfig.case_col ? vocCorrelationConfig.case_col.toString().trim() : 'Case'
    def vocCorrelationSampleTypesRaw = vocCorrelationConfig.sample_types ?: 'Bronchial Brush,Lung Brush'
    def vocCorrelationSampleTypes = vocCorrelationSampleTypesRaw instanceof List ?
        vocCorrelationSampleTypesRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        vocCorrelationSampleTypesRaw.toString().trim()
    boolean vocCorrelationUseLegacySubset = vocCorrelationConfig.containsKey('use_legacy_voc_subset') ? (vocCorrelationConfig.use_legacy_voc_subset as boolean) : true
    def vocCorrelationVocColsRaw = vocCorrelationConfig.voc_columns
    List<String> vocCorrelationVocCols = []
    if( vocCorrelationVocColsRaw instanceof List ) {
        vocCorrelationVocCols = vocCorrelationVocColsRaw.collect { it.toString() }.findAll { it?.trim() }
    } else if( vocCorrelationVocColsRaw ) {
        vocCorrelationVocCols = vocCorrelationVocColsRaw.toString().split(/\r?\n|\|/).collect { it.trim() }.findAll { it }
    }
    def vocCorrelationDirection = vocCorrelationConfig.correlation_direction ? vocCorrelationConfig.correlation_direction.toString().trim().toLowerCase() : 'positive'
    if( !(vocCorrelationDirection in ['positive','negative','both']) ) {
        exit 1, "voc_correlation.correlation_direction must be one of: positive, negative, both"
    }
    def vocCorrelationCasePalette = vocCorrelationConfig.case_palette ?: (indicspeciesGroupPaletteMap[vocCorrelationCaseCol] ?: '')
    def vocCorrelationIsaPalette = vocCorrelationConfig.isa_palette ?: (indicspeciesGroupPaletteMap[vocCorrelationTypeCol] ?: '')

    def measurementAssociationConfig = config.measurement_association ?: [:]
    boolean measurementAssociationRequested = measurementAssociationConfig.containsKey('enabled') ? (measurementAssociationConfig.enabled as boolean) : false
    if( measurementAssociationRequested && !metadataPlotsEnabled ) {
        exit 1, "measurement_association.enabled requires metadata_plots.enabled to be true"
    }
    boolean measurementAssociationEnabled = measurementAssociationRequested
    def measurementAssociationOutputDir = measurementAssociationConfig.output_dir ?: 'measurement_association'
    def measurementAssociationOutputDirAbs = resolveOutputRelative(measurementAssociationOutputDir.toString(), outputDir)
    def measurementAssociationTablePath = measurementAssociationConfig.measurement_table ? resolveOptionalPath(measurementAssociationConfig.measurement_table, configRoot) : null
    if( measurementAssociationEnabled && measurementAssociationTablePath && !new File(measurementAssociationTablePath).exists() ) {
        exit 1, "measurement_association.measurement_table was configured but does not exist: ${measurementAssociationTablePath}"
    }
    def measurementAssociationSampleCol = measurementAssociationConfig.sample_col ? measurementAssociationConfig.sample_col.toString().trim() : metadataPlotsSampleCol
    def measurementAssociationAsvIdCol = measurementAssociationConfig.asv_id_col ? measurementAssociationConfig.asv_id_col.toString().trim() : 'ASV_ID'
    def measurementAssociationMeasurementSampleCol = measurementAssociationConfig.measurement_sample_col ? measurementAssociationConfig.measurement_sample_col.toString().trim() : measurementAssociationSampleCol
    def measurementAssociationMetadataJoinRaw = measurementAssociationConfig.metadata_join_cols ?: ''
    List<String> measurementAssociationMetadataJoinCols = []
    if( measurementAssociationMetadataJoinRaw instanceof List ) {
        measurementAssociationMetadataJoinCols = measurementAssociationMetadataJoinRaw.collect { it.toString().trim() }.findAll { it }
    } else if( measurementAssociationMetadataJoinRaw ) {
        measurementAssociationMetadataJoinCols = measurementAssociationMetadataJoinRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def measurementAssociationMeasurementJoinRaw = measurementAssociationConfig.measurement_join_cols ?: ''
    List<String> measurementAssociationMeasurementJoinCols = []
    if( measurementAssociationMeasurementJoinRaw instanceof List ) {
        measurementAssociationMeasurementJoinCols = measurementAssociationMeasurementJoinRaw.collect { it.toString().trim() }.findAll { it }
    } else if( measurementAssociationMeasurementJoinRaw ) {
        measurementAssociationMeasurementJoinCols = measurementAssociationMeasurementJoinRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def measurementAssociationAliasesRaw = measurementAssociationConfig.measurement_aliases ?: [:]
    if( !(measurementAssociationAliasesRaw instanceof Map) ) {
        exit 1, "measurement_association.measurement_aliases must be a mapping of canonical names to aliases"
    }
    def measurementAssociationAliases = [:]
    measurementAssociationAliasesRaw.each { canonical, aliases ->
        def canonicalName = canonical.toString().trim()
        def aliasValues = aliases instanceof List ? aliases : [aliases]
        def cleanAliases = aliasValues.collect { it.toString().trim() }.findAll { it && it != canonicalName }
        if( canonicalName && cleanAliases ) measurementAssociationAliases[canonicalName] = cleanAliases
    }
    def measurementAssociationAliasesJson = groovy.json.JsonOutput.toJson(measurementAssociationAliases)
    def measurementAssociationColsRaw = measurementAssociationConfig.measurement_cols ?: []
    List<String> measurementAssociationCols = []
    if( measurementAssociationColsRaw instanceof List ) {
        measurementAssociationCols = measurementAssociationColsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( measurementAssociationColsRaw ) {
        measurementAssociationCols = measurementAssociationColsRaw.toString().split(/\r?\n|\|/).collect { it.trim() }.findAll { it }
    }
    def measurementAssociationExcludeRaw = measurementAssociationConfig.exclude_cols ?: []
    List<String> measurementAssociationExcludeCols = []
    if( measurementAssociationExcludeRaw instanceof List ) {
        measurementAssociationExcludeCols = measurementAssociationExcludeRaw.collect { it.toString().trim() }.findAll { it }
    } else if( measurementAssociationExcludeRaw ) {
        measurementAssociationExcludeCols = measurementAssociationExcludeRaw.toString().split(/\r?\n|\|/).collect { it.trim() }.findAll { it }
    }
    def measurementAssociationGroupCol = measurementAssociationConfig.group_col ? measurementAssociationConfig.group_col.toString().trim() : metadataPlotsTypeCol
    def measurementAssociationGroupPalette = measurementAssociationConfig.group_palette ? measurementAssociationConfig.group_palette.toString().trim() : ''
    def measurementAssociationSubsetSource = measurementAssociationConfig.asv_subset_source ?
        measurementAssociationConfig.asv_subset_source.toString().trim().toLowerCase() :
        'internal_filters'
    if( !(measurementAssociationSubsetSource in ['internal_filters', 'spieceasi_standard_filter']) ) {
        exit 1, "measurement_association.asv_subset_source must be one of: internal_filters, spieceasi_standard_filter"
    }
    def measurementAssociationMaxAsvs = measurementAssociationConfig.max_asvs != null ? (measurementAssociationConfig.max_asvs as int) : 0
    def measurementAssociationMinTotal = measurementAssociationConfig.min_total != null ? (measurementAssociationConfig.min_total as double) : 0.0d
    def measurementAssociationMinPrevalence = measurementAssociationConfig.min_prevalence != null ? (measurementAssociationConfig.min_prevalence as double) : 0.0d
    def measurementAssociationTopCorrelations = measurementAssociationConfig.top_correlations ? (measurementAssociationConfig.top_correlations as int) : 100
    def measurementAssociationDirection = measurementAssociationConfig.correlation_direction ? measurementAssociationConfig.correlation_direction.toString().trim().toLowerCase() : 'both'
    if( !(measurementAssociationDirection in ['positive','negative','both']) ) {
        exit 1, "measurement_association.correlation_direction must be one of: positive, negative, both"
    }
    def measurementAssociationMethodsRaw = measurementAssociationConfig.ordination_methods ?: 'cca,rda,dbrda'
    def measurementAssociationMethods = measurementAssociationMethodsRaw instanceof List ?
        measurementAssociationMethodsRaw.collect { it.toString().trim().toLowerCase() }.findAll { it }.join(',') :
        measurementAssociationMethodsRaw.toString().trim().toLowerCase()
    def measurementAssociationOrdinationColsRaw = measurementAssociationConfig.ordination_measurement_cols ?: []
    def measurementAssociationOrdinationCols = measurementAssociationOrdinationColsRaw instanceof List ?
        measurementAssociationOrdinationColsRaw.collect { it.toString().trim() }.findAll { it }.join('|') :
        measurementAssociationOrdinationColsRaw.toString().trim()
    def measurementAssociationPermutations = measurementAssociationConfig.permutations ? (measurementAssociationConfig.permutations as int) : 999
    def measurementAssociationTopVectors = measurementAssociationConfig.top_vectors ? (measurementAssociationConfig.top_vectors as int) : 12
    def measurementAssociationFormatsRaw = measurementAssociationConfig.formats ?: 'pdf,png,svg'
    def measurementAssociationFormats = measurementAssociationFormatsRaw instanceof List ?
        measurementAssociationFormatsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        measurementAssociationFormatsRaw.toString().trim()
    boolean moduleMeasurementAssociationEnabled = measurementAssociationConfig.containsKey('module_analysis_enabled') ?
        (measurementAssociationConfig.module_analysis_enabled as boolean) : false
    def moduleMeasurementSparseRaw = measurementAssociationConfig.sparse_measurement_cols ?: []
    def moduleMeasurementSparseCols = moduleMeasurementSparseRaw instanceof List ?
        moduleMeasurementSparseRaw.collect { it.toString().trim() }.findAll { it }.join('|') :
        moduleMeasurementSparseRaw.toString().trim()
    def moduleMeasurementCruiseCol = measurementAssociationConfig.module_cruise_col ?: 'Cruise'
    def moduleMeasurementDepthCol = measurementAssociationConfig.module_depth_col ?: 'Depth'
    def moduleMeasurementQThreshold = measurementAssociationConfig.module_q_threshold != null ?
        (measurementAssociationConfig.module_q_threshold as double) : 0.05d
    def moduleMeasurementNetworkMetricTopN = measurementAssociationConfig.module_network_metric_top_n != null ?
        (measurementAssociationConfig.module_network_metric_top_n as int) : 10
    if( moduleMeasurementNetworkMetricTopN < 1 ) {
        exit 1, "measurement_association.module_network_metric_top_n must be at least 1"
    }
    def moduleMeasurementPcaScores = measurementAssociationConfig.pca_scores ? resolveOptionalPath(measurementAssociationConfig.pca_scores, configRoot) : null
    def moduleMeasurementPcaLoadings = measurementAssociationConfig.pca_loadings ? resolveOptionalPath(measurementAssociationConfig.pca_loadings, configRoot) : null
    def moduleMeasurementPcaExplained = measurementAssociationConfig.pca_explained ? resolveOptionalPath(measurementAssociationConfig.pca_explained, configRoot) : null
    def moduleMeasurementHybridAssignments = measurementAssociationConfig.hybrid_assignments ? resolveOptionalPath(measurementAssociationConfig.hybrid_assignments, configRoot) : null
    def moduleMeasurementHybridCentroids = measurementAssociationConfig.hybrid_centroids ? resolveOptionalPath(measurementAssociationConfig.hybrid_centroids, configRoot) : null
    def moduleMeasurementPcaInputs = [
        moduleMeasurementPcaScores, moduleMeasurementPcaLoadings,
        moduleMeasurementPcaExplained, moduleMeasurementHybridAssignments,
        moduleMeasurementHybridCentroids,
    ]
    if( moduleMeasurementAssociationEnabled && !measurementAssociationEnabled ) {
        exit 1, "measurement_association.module_analysis_enabled requires measurement_association.enabled=true"
    }
    if( moduleMeasurementPcaInputs.any { it } && !moduleMeasurementPcaInputs.every { it } ) {
        exit 1, "measurement_association module biplot requires pca_scores, pca_loadings, pca_explained, hybrid_assignments, and hybrid_centroids"
    }
    moduleMeasurementPcaInputs.findAll { it }.each { requiredPath ->
        if( !file(requiredPath).exists() ) {
            exit 1, "Module-measurement PCA input not found: ${requiredPath}"
        }
    }

    def titanConfig = config.titan ?: [:]
    boolean titanEnabled = titanConfig.containsKey('enabled') ? (titanConfig.enabled as boolean) : false
    if( titanEnabled && !measurementAssociationEnabled ) {
        exit 1, "titan.enabled requires measurement_association.enabled=true"
    }
    if( titanEnabled && !(config.spieceasi?.enabled as boolean) ) {
        exit 1, "titan.enabled requires spieceasi.enabled=true so the retained_final cohort is available"
    }
    def titanOutputDir = titanConfig.output_dir ?: 'titan'
    def titanOutputDirAbs = resolveOutputRelative(titanOutputDir.toString(), outputDir)
    def titanVariablesRaw = titanConfig.environmental_variables ?: measurementAssociationCols
    def titanVariables = titanVariablesRaw instanceof List ?
        titanVariablesRaw.collect { it.toString().trim() }.findAll { it }.join('|') :
        (titanVariablesRaw ? titanVariablesRaw.toString().trim() : '')
    if( titanEnabled && !titanVariables ) {
        exit 1, "titan.enabled requires titan.environmental_variables or inherited measurement_association.measurement_cols"
    }
    boolean titanTranspose = titanConfig.containsKey('transpose') ?
        (titanConfig.transpose as boolean) :
        (config.spieceasi?.containsKey('transpose') ? (config.spieceasi.transpose as boolean) : true)
    int titanMinSplit = titanConfig.min_split != null ? (titanConfig.min_split as int) : 5
    int titanMinimumSamples = titanConfig.minimum_samples != null ? (titanConfig.minimum_samples as int) : 10
    int titanMinimumOccurrence = titanConfig.minimum_occurrence != null ? (titanConfig.minimum_occurrence as int) : 3
    double titanMinimumPrevalence = titanConfig.minimum_prevalence != null ? (titanConfig.minimum_prevalence as double) : 0.0d
    double titanMinimumMeanRelativeAbundance = titanConfig.minimum_mean_relative_abundance != null ? (titanConfig.minimum_mean_relative_abundance as double) : 0.0d
    int titanPermutations = titanConfig.permutations != null ? (titanConfig.permutations as int) : 250
    int titanBootstrapCount = titanConfig.bootstrap_count != null ? (titanConfig.bootstrap_count as int) : 500
    int titanSeed = titanConfig.seed != null ? (titanConfig.seed as int) : 42
    def titanSourceRepository = (titanConfig.source_repository ?: 'dkahle/TITAN2').toString().trim()
    def titanSourceRevision = (titanConfig.source_revision ?: '0d6c89eda0643c0b1abb74551871adead7a1134d').toString().trim()
    def titanExpectedVersion = (titanConfig.expected_version ?: '2.4.4').toString().trim()
    boolean titanImax = titanConfig.containsKey('imax') ? (titanConfig.imax as boolean) : false
    boolean titanIvTotal = titanConfig.containsKey('iv_total') ? (titanConfig.iv_total as boolean) : false
    double titanPurityCutoff = titanConfig.purity_cutoff != null ? (titanConfig.purity_cutoff as double) : 0.95d
    double titanReliabilityCutoff = titanConfig.reliability_cutoff != null ? (titanConfig.reliability_cutoff as double) : 0.95d
    int titanNcpus = titanConfig.ncpus != null ? (titanConfig.ncpus as int) : pipelineThreads
    boolean titanMemory = titanConfig.containsKey('memory') ? (titanConfig.memory as boolean) : false
    boolean titanPlotEnabled = titanConfig.containsKey('plot_enabled') ? (titanConfig.plot_enabled as boolean) : true
    int titanRankedTopN = titanConfig.ranked_plot_top_n != null ? (titanConfig.ranked_plot_top_n as int) : 50
    def titanFormatsRaw = titanConfig.formats ?: measurementAssociationFormats
    def titanFormats = titanFormatsRaw instanceof List ?
        titanFormatsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        titanFormatsRaw.toString().trim()
    if( titanMinSplit < 3 ) exit 1, "titan.min_split must be at least 3"
    if( titanMinimumSamples < 10 ) exit 1, "titan.minimum_samples cannot be below TITAN2's hard minimum of 10"
    if( titanMinimumOccurrence < 3 ) exit 1, "titan.minimum_occurrence cannot be below TITAN2's hard minimum of 3"
    if( titanMinimumPrevalence < 0 || titanMinimumPrevalence > 1 || titanMinimumMeanRelativeAbundance < 0 || titanMinimumMeanRelativeAbundance > 1 ) {
        exit 1, "titan minimum_prevalence and minimum_mean_relative_abundance must be between 0 and 1"
    }
    if( titanPermutations < 1 || titanBootstrapCount < 1 || titanNcpus < 1 ) {
        exit 1, "titan permutations, bootstrap_count, and ncpus must be positive"
    }
    if( titanNcpus > pipelineThreads ) {
        exit 1, "titan.ncpus (${titanNcpus}) cannot exceed resources.threads (${pipelineThreads})"
    }
    if( titanPurityCutoff < 0 || titanPurityCutoff > 1 || titanReliabilityCutoff < 0 || titanReliabilityCutoff > 1 ) {
        exit 1, "titan purity_cutoff and reliability_cutoff must be between 0 and 1"
    }
    if( titanRankedTopN < 0 ) exit 1, "titan.ranked_plot_top_n must be nonnegative; 0 plots all taxa"
    if( !(titanSourceRepository ==~ /[A-Za-z0-9_.-]+\/[A-Za-z0-9_.-]+/) ) {
        exit 1, "titan.source_repository must use GitHub owner/repository syntax"
    }
    if( !(titanSourceRevision ==~ /[0-9a-fA-F]{40}/) ) {
        exit 1, "titan.source_revision must be a full 40-character Git commit SHA"
    }
    if( !titanExpectedVersion ) exit 1, "titan.expected_version cannot be empty"

    def microbialCompartmentConfig = config.microbial_compartments ?: [:]
    boolean microbialCompartmentEnabled = microbialCompartmentConfig.containsKey('enabled') ? (microbialCompartmentConfig.enabled as boolean) : false
    if( microbialCompartmentEnabled && !metadataPlotsEnabled ) {
        exit 1, "microbial_compartments.enabled requires metadata_plots.enabled=true"
    }
    if( microbialCompartmentEnabled && !(config.spieceasi?.enabled as boolean) ) {
        exit 1, "microbial_compartments.enabled requires spieceasi.enabled=true so retained_final is available"
    }
    def microbialCompartmentOutputDir = microbialCompartmentConfig.output_dir ?: 'microbial_compartments'
    def microbialCompartmentOutputDirAbs = resolveOutputRelative(microbialCompartmentOutputDir.toString(), outputDir)
    def microbialCompartmentSampleCol = microbialCompartmentConfig.sample_col ?: metadataPlotsSampleCol
    boolean microbialCompartmentTranspose = microbialCompartmentConfig.containsKey('transpose') ?
        (microbialCompartmentConfig.transpose as boolean) :
        (config.spieceasi?.containsKey('transpose') ? (config.spieceasi.transpose as boolean) : true)
    double microbialCompartmentMinimumPrevalence = microbialCompartmentConfig.minimum_prevalence != null ? (microbialCompartmentConfig.minimum_prevalence as double) : 0.0d
    double microbialCompartmentMinimumMaxRelativeAbundance = microbialCompartmentConfig.minimum_max_relative_abundance != null ? (microbialCompartmentConfig.minimum_max_relative_abundance as double) : 0.0d
    def microbialCompartmentZeroReplacement = microbialCompartmentConfig.zero_replacement ? microbialCompartmentConfig.zero_replacement.toString().trim().toLowerCase() : 'multiplicative'
    double microbialCompartmentMultiplicativeDelta = microbialCompartmentConfig.multiplicative_delta != null ? (microbialCompartmentConfig.multiplicative_delta as double) : 0.0d
    double microbialCompartmentPseudocount = microbialCompartmentConfig.pseudocount != null ? (microbialCompartmentConfig.pseudocount as double) : 1e-6d
    int microbialCompartmentKMin = microbialCompartmentConfig.k_min != null ? (microbialCompartmentConfig.k_min as int) : 2
    int microbialCompartmentKMax = microbialCompartmentConfig.k_max != null ? (microbialCompartmentConfig.k_max as int) : 10
    int microbialCompartmentMinClusterSize = microbialCompartmentConfig.min_cluster_size != null ? (microbialCompartmentConfig.min_cluster_size as int) : 5
    double microbialCompartmentMinClusterFraction = microbialCompartmentConfig.min_cluster_fraction != null ? (microbialCompartmentConfig.min_cluster_fraction as double) : 0.02d
    double microbialCompartmentMinMeanSilhouette = microbialCompartmentConfig.min_mean_silhouette != null ? (microbialCompartmentConfig.min_mean_silhouette as double) : 0.25d
    double microbialCompartmentMinStabilityAri = microbialCompartmentConfig.min_stability_ari != null ? (microbialCompartmentConfig.min_stability_ari as double) : 0.70d
    double microbialCompartmentMinClusterJaccard = microbialCompartmentConfig.min_cluster_jaccard != null ? (microbialCompartmentConfig.min_cluster_jaccard as double) : 0.75d
    int microbialCompartmentStabilityReplicates = microbialCompartmentConfig.stability_replicates != null ? (microbialCompartmentConfig.stability_replicates as int) : 200
    double microbialCompartmentStabilitySampleFraction = microbialCompartmentConfig.stability_sample_fraction != null ? (microbialCompartmentConfig.stability_sample_fraction as double) : 0.80d
    def microbialCompartmentStabilityBlockCol = (microbialCompartmentConfig.stability_block_col ?: '').toString().trim()
    def microbialCompartmentStabilityStratumCol = (microbialCompartmentConfig.stability_stratum_col ?: '').toString().trim()
    double microbialCompartmentStabilityPrimaryQuantile = microbialCompartmentConfig.stability_primary_quantile != null ? (microbialCompartmentConfig.stability_primary_quantile as double) : 0.25d
    double microbialCompartmentStabilityNearTieTolerance = microbialCompartmentConfig.stability_near_tie_tolerance != null ? (microbialCompartmentConfig.stability_near_tie_tolerance as double) : 0.02d
    boolean microbialCompartmentBlockedRobustnessEnabled = microbialCompartmentConfig.containsKey('blocked_robustness_enabled') ? (microbialCompartmentConfig.blocked_robustness_enabled as boolean) : false
    boolean microbialCompartmentPredictionStrengthEnabled = microbialCompartmentConfig.containsKey('prediction_strength_enabled') ? (microbialCompartmentConfig.prediction_strength_enabled as boolean) : true
    boolean microbialCompartmentHierarchicalEnabled = microbialCompartmentConfig.containsKey('hierarchical_enabled') ? (microbialCompartmentConfig.hierarchical_enabled as boolean) : true
    int microbialCompartmentSeed = microbialCompartmentConfig.seed != null ? (microbialCompartmentConfig.seed as int) : 42
    int microbialCompartmentNcpus = microbialCompartmentConfig.ncpus != null ? (microbialCompartmentConfig.ncpus as int) : pipelineThreads
    def microbialCompartmentEnvironmentalColsRaw = microbialCompartmentConfig.environmental_compartment_cols ?: []
    def microbialCompartmentEnvironmentalCols = microbialCompartmentEnvironmentalColsRaw instanceof List ?
        microbialCompartmentEnvironmentalColsRaw.collect { it.toString().trim() }.findAll { it }.join('|') :
        microbialCompartmentEnvironmentalColsRaw.toString().trim()
    def microbialCompartmentDepthCol = microbialCompartmentConfig.depth_col ?: 'Depth'
    def microbialCompartmentDateCol = microbialCompartmentConfig.date_col ?: 'date'
    int microbialCompartmentDominantTopN = microbialCompartmentConfig.dominant_top_n != null ? (microbialCompartmentConfig.dominant_top_n as int) : 10
    boolean microbialCompartmentPlotEnabled = microbialCompartmentConfig.containsKey('plot_enabled') ? (microbialCompartmentConfig.plot_enabled as boolean) : true
    def microbialCompartmentFormatsRaw = microbialCompartmentConfig.formats ?: 'pdf,png,svg'
    def microbialCompartmentFormats = microbialCompartmentFormatsRaw instanceof List ?
        microbialCompartmentFormatsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        microbialCompartmentFormatsRaw.toString().trim()
    if( microbialCompartmentMinimumPrevalence < 0 || microbialCompartmentMinimumPrevalence > 1 ||
        microbialCompartmentMinimumMaxRelativeAbundance < 0 || microbialCompartmentMinimumMaxRelativeAbundance > 1 ) {
        exit 1, "microbial_compartments filtering thresholds must be between 0 and 1"
    }
    if( !(microbialCompartmentZeroReplacement in ['multiplicative', 'pseudocount']) ) {
        exit 1, "microbial_compartments.zero_replacement must be multiplicative or pseudocount"
    }
    if( microbialCompartmentMultiplicativeDelta < 0 || microbialCompartmentPseudocount <= 0 ) {
        exit 1, "microbial_compartments multiplicative_delta must be nonnegative and pseudocount must be positive"
    }
    if( microbialCompartmentKMin < 2 || microbialCompartmentKMax < microbialCompartmentKMin ) {
        exit 1, "microbial_compartments requires 2 <= k_min <= k_max"
    }
    if( microbialCompartmentMinClusterSize < 2 || microbialCompartmentMinClusterFraction < 0 || microbialCompartmentMinClusterFraction >= 1 ) {
        exit 1, "microbial_compartments cluster-size constraints are invalid"
    }
    if( microbialCompartmentMinMeanSilhouette < -1 || microbialCompartmentMinMeanSilhouette > 1 ||
        microbialCompartmentMinStabilityAri < -1 || microbialCompartmentMinStabilityAri > 1 ||
        microbialCompartmentMinClusterJaccard < 0 || microbialCompartmentMinClusterJaccard > 1 ) {
        exit 1, "microbial_compartments silhouette, ARI, and Jaccard thresholds are invalid"
    }
    if( microbialCompartmentStabilityPrimaryQuantile <= 0 || microbialCompartmentStabilityPrimaryQuantile >= 0.5 ||
        microbialCompartmentStabilityNearTieTolerance < 0 || microbialCompartmentStabilityNearTieTolerance > 1 ) {
        exit 1, "microbial_compartments stability quantile or near-tie tolerance is invalid"
    }
    if( !microbialCompartmentStabilityBlockCol ) {
        exit 1, "microbial_compartments.stability_block_col is required for within-block stability"
    }
    if( microbialCompartmentBlockedRobustnessEnabled && !microbialCompartmentStabilityStratumCol ) {
        exit 1, "microbial_compartments.stability_stratum_col is required when blocked robustness is enabled"
    }
    boolean indicspeciesRequiresMicrobialCompartments = indicspeciesGroupCols.contains('microbial_compartment') ||
        indicspeciesBlockedCols.split(',').collect { it.trim() }.contains('microbial_compartment') ||
        indicspeciesStratifiedSpecs.any { spec -> spec.split(/::/, -1).take(2).contains('microbial_compartment') }
    boolean indicspeciesRunEarly = indicspeciesEnabled && !indicspeciesRequiresMicrobialCompartments
    if( indicspeciesRequiresMicrobialCompartments && !microbialCompartmentEnabled ) {
        exit 1, "Indicator-species analyses using microbial_compartment require microbial_compartments.enabled=true"
    }
    if( microbialCompartmentStabilityReplicates < 1 || microbialCompartmentStabilitySampleFraction <= 0 ||
        microbialCompartmentStabilitySampleFraction > 1 || microbialCompartmentNcpus < 1 ||
        microbialCompartmentNcpus > pipelineThreads || microbialCompartmentDominantTopN < 1 ) {
        exit 1, "microbial_compartments replicate, sampling, CPU, or dominant-ASV settings are invalid"
    }

    def microbialStateInterpretationConfig = config.microbial_state_interpretation ?: [:]
    boolean microbialStateInterpretationEnabled = microbialStateInterpretationConfig.containsKey('enabled') ? (microbialStateInterpretationConfig.enabled as boolean) : false
    if( microbialStateInterpretationEnabled && (!microbialCompartmentEnabled || !(config.spieceasi?.modules_enabled as boolean)) ) {
        exit 1, "microbial_state_interpretation.enabled requires microbial_compartments.enabled and spieceasi.modules_enabled"
    }
    def microbialStateInterpretationOutputDir = microbialStateInterpretationConfig.output_dir ?: 'microbial_state_interpretation'
    def microbialStateInterpretationOutputDirAbs = resolveOutputRelative(microbialStateInterpretationOutputDir.toString(), outputDir)
    def microbialStateInterpretationSampleCol = microbialStateInterpretationConfig.sample_col ?: microbialCompartmentSampleCol
    def microbialStateInterpretationCruiseCol = microbialStateInterpretationConfig.cruise_col ?: 'Cruise'
    def microbialStateInterpretationDepthCol = microbialStateInterpretationConfig.depth_col ?: microbialCompartmentDepthCol
    def microbialStateInterpretationSeasonCol = microbialStateInterpretationConfig.season_col ?: 'Season'
    def microbialStateInterpretationYearCol = microbialStateInterpretationConfig.year_col ?: 'Year'
    def microbialStateInterpretationDateCol = microbialStateInterpretationConfig.date_col ?: microbialCompartmentDateCol
    def microbialStateInterpretationRenewalCol = microbialStateInterpretationConfig.renewal_col ?: 'renewal_phase'
    def microbialStateInterpretationRenewalLevelsRaw = microbialStateInterpretationConfig.renewal_levels ?: ['baseline','renewal','post-renewal']
    def microbialStateInterpretationRenewalLevels = microbialStateInterpretationRenewalLevelsRaw instanceof List ? microbialStateInterpretationRenewalLevelsRaw.join(',') : microbialStateInterpretationRenewalLevelsRaw.toString()
    def microbialStateInterpretationEnvironmentalColsRaw = microbialStateInterpretationConfig.environmental_compartment_cols ?: microbialCompartmentEnvironmentalColsRaw
    def microbialStateInterpretationEnvironmentalCols = microbialStateInterpretationEnvironmentalColsRaw instanceof List ?
        microbialStateInterpretationEnvironmentalColsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        microbialStateInterpretationEnvironmentalColsRaw.toString().trim()
    def microbialStateInterpretationHybridCol = microbialStateInterpretationConfig.hybrid_col ?: 'o2_subcompartment_final'
    int microbialStateInterpretationPermutations = microbialStateInterpretationConfig.permutations != null ? (microbialStateInterpretationConfig.permutations as int) : 999
    int microbialStateInterpretationBootstrapReplicates = microbialStateInterpretationConfig.bootstrap_replicates != null ? (microbialStateInterpretationConfig.bootstrap_replicates as int) : 5000
    int microbialStateInterpretationSeed = microbialStateInterpretationConfig.seed != null ? (microbialStateInterpretationConfig.seed as int) : 42
    def microbialStateInterpretationHybridPalette = microbialStateInterpretationConfig.hybrid_palette ?: ''
    def microbialStateInterpretationHybridOrderRaw = microbialStateInterpretationConfig.hybrid_order ?: []
    def microbialStateInterpretationHybridOrder = microbialStateInterpretationHybridOrderRaw instanceof List ? microbialStateInterpretationHybridOrderRaw.join(',') : microbialStateInterpretationHybridOrderRaw.toString()
    def microbialStateInterpretationMcPalette = microbialStateInterpretationConfig.mc_palette ?: 'MC1=#0072B2,MC2=#E69F00,MC3=#009E73,MC4=#CC79A7'
    def microbialStateInterpretationDepthPalette = microbialStateInterpretationConfig.depth_palette ?: ''
    def microbialStateInterpretationMcOrderRaw = microbialStateInterpretationConfig.mc_order ?: ['MC1','MC2','MC3','MC4']
    def microbialStateInterpretationMcOrder = microbialStateInterpretationMcOrderRaw instanceof List ? microbialStateInterpretationMcOrderRaw.join(',') : microbialStateInterpretationMcOrderRaw.toString()
    int microbialStateInterpretationAsvTopN = microbialStateInterpretationConfig.asv_top_n != null ? (microbialStateInterpretationConfig.asv_top_n as int) : 20
    def microbialStateInterpretationSelectedAsvsRaw = microbialStateInterpretationConfig.selected_asvs ?: []
    def microbialStateInterpretationSelectedAsvs = microbialStateInterpretationSelectedAsvsRaw instanceof List ? microbialStateInterpretationSelectedAsvsRaw.join(',') : microbialStateInterpretationSelectedAsvsRaw.toString()
    def microbialStateInterpretationAsvSelectionMetric = microbialStateInterpretationConfig.asv_selection_metric ?: 'mean_relative_abundance'
    double microbialStateInterpretationMaximumDepth = microbialStateInterpretationConfig.maximum_depth != null ? (microbialStateInterpretationConfig.maximum_depth as double) : 210.0d
    double microbialStateInterpretationDepthStep = microbialStateInterpretationConfig.depth_step != null ? (microbialStateInterpretationConfig.depth_step as double) : 1.0d
    int microbialStateInterpretationTimeSubdivisions = microbialStateInterpretationConfig.time_subdivisions_per_month != null ? (microbialStateInterpretationConfig.time_subdivisions_per_month as int) : 4
    double microbialStateInterpretationTimeSigma = microbialStateInterpretationConfig.time_sigma_months != null ? (microbialStateInterpretationConfig.time_sigma_months as double) : 0.75d
    double microbialStateInterpretationDepthSigma = microbialStateInterpretationConfig.depth_visual_sigma_m != null ? (microbialStateInterpretationConfig.depth_visual_sigma_m as double) : 5.0d
    int microbialStateInterpretationContourTimeStepDays = microbialStateInterpretationConfig.contour_time_step_days != null ? (microbialStateInterpretationConfig.contour_time_step_days as int) : 7
    double microbialStateInterpretationContourMaxTimeSupportDays = microbialStateInterpretationConfig.contour_max_time_support_days != null ? (microbialStateInterpretationConfig.contour_max_time_support_days as double) : 90.0d
    double microbialStateInterpretationContourMaxDepthSupportM = microbialStateInterpretationConfig.contour_max_depth_support_m != null ? (microbialStateInterpretationConfig.contour_max_depth_support_m as double) : 30.0d
    def microbialStateInterpretationExcludeCurtainPattern = microbialStateInterpretationConfig.exclude_curtain_label_pattern ?: '(?i)(?:other|outlier)'
    def microbialStateInterpretationRenewalEvents = microbialStateInterpretationConfig.renewal_events ?
        resolveOptionalPath(microbialStateInterpretationConfig.renewal_events.toString(), configRoot) : null
    def microbialStateInterpretationRenewalDateCol = (microbialStateInterpretationConfig.renewal_date_col ?: 'start_date').toString()
    def microbialStateInterpretationFormatsRaw = microbialStateInterpretationConfig.formats ?: 'pdf,png,svg'
    def microbialStateInterpretationFormats = microbialStateInterpretationFormatsRaw instanceof List ? microbialStateInterpretationFormatsRaw.join(',') : microbialStateInterpretationFormatsRaw.toString()
    if( microbialStateInterpretationEnabled && !microbialStateInterpretationHybridPalette ) exit 1, "microbial_state_interpretation.hybrid_palette is required"
    if( microbialStateInterpretationPermutations < 1 || microbialStateInterpretationBootstrapReplicates < 1 || microbialStateInterpretationAsvTopN < 1 ||
        microbialStateInterpretationMaximumDepth <= 0 || microbialStateInterpretationDepthStep <= 0 || microbialStateInterpretationTimeSubdivisions < 1 ||
        microbialStateInterpretationContourTimeStepDays < 1 || microbialStateInterpretationContourMaxTimeSupportDays <= 0 || microbialStateInterpretationContourMaxDepthSupportM <= 0 ) {
        exit 1, "microbial_state_interpretation replicate, ASV, or curtain settings are invalid"
    }

    def ecologicalContextAtlasConfig = config.ecological_context_atlas ?: [:]
    boolean ecologicalContextAtlasEnabled = ecologicalContextAtlasConfig.containsKey('enabled') ? (ecologicalContextAtlasConfig.enabled as boolean) : false
    if( ecologicalContextAtlasEnabled && (!titanEnabled || !measurementAssociationEnabled || !microbialCompartmentEnabled || !(config.spieceasi?.modules_enabled as boolean) || !(config.spieceasi?.network_enabled as boolean) || !(config.asv_mag_network?.enabled as boolean) || !(config.group_guild_function?.enabled as boolean)) ) {
        exit 1, "ecological_context_atlas.enabled requires TITAN, measurement association, microbial compartments, SPIEC-EASI modules, network analysis, ASV-MAG network, and group-guild-function module associations"
    }
    def ecologicalContextAtlasOutputDir = ecologicalContextAtlasConfig.output_dir ?: 'ecological_context_atlas'
    def ecologicalContextAtlasOutputDirAbs = resolveOutputRelative(ecologicalContextAtlasOutputDir.toString(), outputDir)
    def ecologicalContextAtlasSampleCol = ecologicalContextAtlasConfig.sample_col ?: microbialCompartmentSampleCol
    def ecologicalContextAtlasCruiseCol = ecologicalContextAtlasConfig.cruise_col ?: microbialStateInterpretationCruiseCol
    def ecologicalContextAtlasDepthCol = ecologicalContextAtlasConfig.depth_col ?: microbialCompartmentDepthCol
    def ecologicalContextAtlasDateCol = ecologicalContextAtlasConfig.date_col ?: microbialCompartmentDateCol
    def ecologicalContextAtlasContextColsRaw = ecologicalContextAtlasConfig.context_cols ?: ['o2_compartment','gmm_component','o2_subcompartment_final','microbial_compartment','Season','renewal_phase','cruise_group']
    def ecologicalContextAtlasContextCols = ecologicalContextAtlasContextColsRaw instanceof List ? ecologicalContextAtlasContextColsRaw.join(',') : ecologicalContextAtlasContextColsRaw.toString()
    // Empty or omitted means every TITAN-configured measurement. Variables
    // skipped by TITAN for absent/insufficient observations are omitted by the
    // atlas plotter rather than rendered as empty networks.
    def ecologicalContextAtlasNetworkVariablesRaw = ecologicalContextAtlasConfig.network_variables ?: titanVariablesRaw
    def ecologicalContextAtlasNetworkVariables = ecologicalContextAtlasNetworkVariablesRaw instanceof List ? ecologicalContextAtlasNetworkVariablesRaw.join(',') : ecologicalContextAtlasNetworkVariablesRaw.toString()
    def ecologicalContextAtlasLinkageQThreshold = ecologicalContextAtlasConfig.linkage_q_threshold != null ? (ecologicalContextAtlasConfig.linkage_q_threshold as double) : 0.05d
    def ecologicalContextAtlasSelectedAsvsRaw = ecologicalContextAtlasConfig.selected_asvs ?: []
    def ecologicalContextAtlasSelectedAsvs = ecologicalContextAtlasSelectedAsvsRaw instanceof List ? ecologicalContextAtlasSelectedAsvsRaw.join(',') : ecologicalContextAtlasSelectedAsvsRaw.toString()
    def ecologicalContextAtlasContextSheetSource = ecologicalContextAtlasConfig.context_sheet_source ?: 'mag_linked'
    def ecologicalContextAtlasRnaDnaLog2TpmRatio = ecologicalContextAtlasConfig.rna_dna_log2_tpm_ratio ? resolveOptionalPath(ecologicalContextAtlasConfig.rna_dna_log2_tpm_ratio, configRoot) : null
    def ecologicalContextAtlasRnaDnaLog2TpmSe = ecologicalContextAtlasConfig.rna_dna_log2_tpm_se ? resolveOptionalPath(ecologicalContextAtlasConfig.rna_dna_log2_tpm_se, configRoot) : null
    def ecologicalContextAtlasRnaDnaMaxLog2Se = ecologicalContextAtlasConfig.rna_dna_max_log2_se != null ? (ecologicalContextAtlasConfig.rna_dna_max_log2_se as double) : 0.3d
    def ecologicalContextAtlasRnaDnaMagIdMode = ecologicalContextAtlasConfig.rna_dna_mag_id_mode ? ecologicalContextAtlasConfig.rna_dna_mag_id_mode.toString().trim().toLowerCase() : (config.asv_mag_network?.mag_id_mode ?: 'exact').toString().trim().toLowerCase()
    def ecologicalContextAtlasFormatsRaw = ecologicalContextAtlasConfig.formats ?: 'pdf,png,svg'
    def ecologicalContextAtlasFormats = ecologicalContextAtlasFormatsRaw instanceof List ? ecologicalContextAtlasFormatsRaw.join(',') : ecologicalContextAtlasFormatsRaw.toString()
    if( !(ecologicalContextAtlasContextSheetSource in ['explicit','mag_linked']) ) exit 1, "ecological_context_atlas.context_sheet_source must be explicit or mag_linked"
    if( ecologicalContextAtlasRnaDnaMaxLog2Se < 0 ) exit 1, "ecological_context_atlas.rna_dna_max_log2_se must be nonnegative"
    if( ecologicalContextAtlasLinkageQThreshold <= 0 || ecologicalContextAtlasLinkageQThreshold > 1 ) exit 1, "ecological_context_atlas.linkage_q_threshold must be in (0, 1]"
    if( !(ecologicalContextAtlasRnaDnaMagIdMode in ['exact','suffix_after_double_underscore']) ) exit 1, "ecological_context_atlas.rna_dna_mag_id_mode must be exact or suffix_after_double_underscore"
    ['rna_dna_log2_tpm_ratio': ecologicalContextAtlasRnaDnaLog2TpmRatio, 'rna_dna_log2_tpm_se': ecologicalContextAtlasRnaDnaLog2TpmSe].each { label, candidate ->
        if( candidate && !new File(candidate).isFile() ) exit 1, "ecological_context_atlas.${label} not found: ${candidate}"
    }

    def communityTurnoverConfig = config.community_turnover ?: [:]
    boolean communityTurnoverEnabled = communityTurnoverConfig.containsKey('enabled') ? (communityTurnoverConfig.enabled as boolean) : false
    if( communityTurnoverEnabled && !metadataPlotsEnabled ) {
        exit 1, "community_turnover.enabled requires metadata_plots.enabled=true"
    }
    if( communityTurnoverEnabled && !(config.spieceasi?.enabled as boolean) ) {
        exit 1, "community_turnover.enabled requires spieceasi.enabled=true so retained_final is available"
    }
    def communityTurnoverOutputDir = communityTurnoverConfig.output_dir ?: 'community_turnover'
    def communityTurnoverOutputDirAbs = resolveOutputRelative(communityTurnoverOutputDir.toString(), outputDir)
    def communityTurnoverSampleCol = communityTurnoverConfig.sample_col ?: metadataPlotsSampleCol
    def communityTurnoverProfileCol = communityTurnoverConfig.profile_col ?: 'Cruise'
    def communityTurnoverDateCol = communityTurnoverConfig.date_col ?: 'date'
    def communityTurnoverDepthCol = communityTurnoverConfig.depth_col ?: 'Depth'
    boolean communityTurnoverTranspose = communityTurnoverConfig.containsKey('transpose') ?
        (communityTurnoverConfig.transpose as boolean) :
        (config.spieceasi?.containsKey('transpose') ? (config.spieceasi.transpose as boolean) : true)
    double communityTurnoverMinimumPrevalence = communityTurnoverConfig.minimum_prevalence != null ? (communityTurnoverConfig.minimum_prevalence as double) : 0.0d
    double communityTurnoverMinimumMaxRelativeAbundance = communityTurnoverConfig.minimum_max_relative_abundance != null ? (communityTurnoverConfig.minimum_max_relative_abundance as double) : 0.0d
    def communityTurnoverZeroReplacement = communityTurnoverConfig.zero_replacement ? communityTurnoverConfig.zero_replacement.toString().trim().toLowerCase() : 'multiplicative'
    double communityTurnoverMultiplicativeDelta = communityTurnoverConfig.multiplicative_delta != null ? (communityTurnoverConfig.multiplicative_delta as double) : 0.0d
    double communityTurnoverPseudocount = communityTurnoverConfig.pseudocount != null ? (communityTurnoverConfig.pseudocount as double) : 1e-6d
    def communityTurnoverDistanceMetricsRaw = communityTurnoverConfig.distance_metrics ?: ['aitchison', 'braycurtis']
    def communityTurnoverDistanceMetrics = communityTurnoverDistanceMetricsRaw instanceof List ?
        communityTurnoverDistanceMetricsRaw.collect { it.toString().trim().toLowerCase() }.findAll { it }.join(',') :
        communityTurnoverDistanceMetricsRaw.toString().trim().toLowerCase()
    def communityTurnoverPrimaryMetric = communityTurnoverConfig.primary_metric ? communityTurnoverConfig.primary_metric.toString().trim().toLowerCase() : 'aitchison'
    double communityTurnoverFixedDepthMinProfileFraction = communityTurnoverConfig.fixed_depth_min_profile_fraction != null ? (communityTurnoverConfig.fixed_depth_min_profile_fraction as double) : 0.50d
    double communityTurnoverMinimumTimeDifferenceDays = communityTurnoverConfig.minimum_time_difference_days != null ? (communityTurnoverConfig.minimum_time_difference_days as double) : 0.0d
    def communityTurnoverEnvironmentalColsRaw = communityTurnoverConfig.environmental_compartment_cols ?: []
    def communityTurnoverEnvironmentalCols = communityTurnoverEnvironmentalColsRaw instanceof List ?
        communityTurnoverEnvironmentalColsRaw.collect { it.toString().trim() }.findAll { it }.join('|') :
        communityTurnoverEnvironmentalColsRaw.toString().trim()
    def communityTurnoverPrimaryEnvironmentalCol = communityTurnoverConfig.primary_environmental_compartment_col ?: ''
    int communityTurnoverLcbdPermutations = communityTurnoverConfig.lcbd_permutations != null ? (communityTurnoverConfig.lcbd_permutations as int) : 999
    int communityTurnoverBoundaryPermutations = communityTurnoverConfig.boundary_permutations != null ? (communityTurnoverConfig.boundary_permutations as int) : 9999
    int communityTurnoverBoundaryBootstrapReplicates = communityTurnoverConfig.boundary_bootstrap_replicates != null ? (communityTurnoverConfig.boundary_bootstrap_replicates as int) : 5000
    boolean communityTurnoverWithinProfileLcbdEnabled = communityTurnoverConfig.containsKey('within_profile_lcbd_enabled') ? (communityTurnoverConfig.within_profile_lcbd_enabled as boolean) : true
    int communityTurnoverWithinProfileLcbdMinSamples = communityTurnoverConfig.within_profile_lcbd_min_samples != null ? (communityTurnoverConfig.within_profile_lcbd_min_samples as int) : 3
    int communityTurnoverSeed = communityTurnoverConfig.seed != null ? (communityTurnoverConfig.seed as int) : 42
    boolean communityTurnoverPlotEnabled = communityTurnoverConfig.containsKey('plot_enabled') ? (communityTurnoverConfig.plot_enabled as boolean) : true
    def communityTurnoverFormatsRaw = communityTurnoverConfig.formats ?: 'pdf,png,svg'
    def communityTurnoverFormats = communityTurnoverFormatsRaw instanceof List ?
        communityTurnoverFormatsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        communityTurnoverFormatsRaw.toString().trim()
    def resolvedTurnoverMetrics = communityTurnoverDistanceMetrics.split(',').collect { it.trim() }.findAll { it }
    if( resolvedTurnoverMetrics.isEmpty() || resolvedTurnoverMetrics.any { !(it in ['aitchison', 'braycurtis']) } || !(communityTurnoverPrimaryMetric in resolvedTurnoverMetrics) ) {
        exit 1, "community_turnover distance_metrics must contain aitchison and/or braycurtis and include primary_metric"
    }
    if( communityTurnoverMinimumPrevalence < 0 || communityTurnoverMinimumPrevalence > 1 ||
        communityTurnoverMinimumMaxRelativeAbundance < 0 || communityTurnoverMinimumMaxRelativeAbundance > 1 ||
        communityTurnoverFixedDepthMinProfileFraction <= 0 || communityTurnoverFixedDepthMinProfileFraction > 1 ) {
        exit 1, "community_turnover prevalence, abundance, and common-depth fractions must be within valid bounds"
    }
    if( !(communityTurnoverZeroReplacement in ['multiplicative', 'pseudocount']) ||
        communityTurnoverMultiplicativeDelta < 0 || communityTurnoverPseudocount <= 0 ||
        communityTurnoverMinimumTimeDifferenceDays < 0 || communityTurnoverLcbdPermutations < 1 ||
        communityTurnoverBoundaryPermutations < 1 || communityTurnoverBoundaryBootstrapReplicates < 100 ||
        communityTurnoverWithinProfileLcbdMinSamples < 3 ) {
        exit 1, "community_turnover replacement, temporal, or LCBD settings are invalid"
    }

    def communityPredictorConfig = config.community_predictor_comparison ?: [:]
    boolean communityPredictorRequested = communityPredictorConfig.containsKey('enabled') ? (communityPredictorConfig.enabled as boolean) : false
    if( communityPredictorRequested && !metadataPlotsEnabled ) {
        exit 1, "community_predictor_comparison.enabled requires metadata_plots.enabled to be true"
    }
    boolean communityPredictorEnabled = communityPredictorRequested
    def communityPredictorOutputDir = communityPredictorConfig.output_dir ?: 'community_predictor_comparison'
    def communityPredictorOutputDirAbs = resolveOutputRelative(communityPredictorOutputDir.toString(), outputDir)
    def communityPredictorSampleCol = communityPredictorConfig.sample_col ?: metadataPlotsSampleCol
    def communityPredictorCruiseCol = communityPredictorConfig.cruise_col ?: 'Cruise'
    def communityPredictorYearCol = communityPredictorConfig.year_col ?: 'Year'
    def communityPredictorDateCol = communityPredictorConfig.date_col ?: 'date'
    def communityPredictorSeasonCol = communityPredictorConfig.season_col ?: 'Season'
    def communityPredictorDepthCol = communityPredictorConfig.depth_col ?: 'Depth'
    double communityPredictorCruiseDepthMinPrevalence = communityPredictorConfig.cruise_depth_min_prevalence != null ? (communityPredictorConfig.cruise_depth_min_prevalence as double) : 0.50
    def communityPredictorPeaCol = communityPredictorConfig.pea_col ?: 'pea_J_m3'
    def communityPredictorCentroidCol = communityPredictorConfig.centroid_col ?: 'depth_centroid_distance'
    def communityPredictorCruiseGroupCol = communityPredictorConfig.cruise_group_col ?: 'cruise_group'
    def communityPredictorCruiseGroupProbabilityCol = communityPredictorConfig.cruise_group_probability_col ?: 'cruise_group_max_prob'
    def communityPredictorCruiseGroupUncertainCol = communityPredictorConfig.cruise_group_uncertain_col ?: 'cruise_group_uncertain'
    def communityPredictorRenewalGroupCol = communityPredictorConfig.renewal_group_col ?: 'renewal_phase'
    def communityPredictorO2Col = communityPredictorConfig.o2_group_col ?: 'o2_compartment'
    def communityPredictorGmmCol = communityPredictorConfig.gmm_group_col ?: 'gmm_component'
    def communityPredictorHybridCol = communityPredictorConfig.hybrid_group_col ?: 'o2_subcompartment_final'
    def communityPredictorSourceCol = communityPredictorConfig.assignment_source_col ?: 'o2_subcompartment_final_assignment_source'
    def communityPredictorObservedLabel = communityPredictorConfig.observed_source_label ?: 'observed'
    def communityPredictorMinGroupN = communityPredictorConfig.min_group_n ? (communityPredictorConfig.min_group_n as int) : 3
    def communityPredictorPermutations = communityPredictorConfig.permutations ? (communityPredictorConfig.permutations as int) : 999
    def communityPredictorSeed = communityPredictorConfig.seed ? (communityPredictorConfig.seed as int) : 42
    def communityPredictorFormatsRaw = communityPredictorConfig.formats ?: 'pdf,png,svg'
    def communityPredictorFormats = communityPredictorFormatsRaw instanceof List ?
        communityPredictorFormatsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        communityPredictorFormatsRaw.toString().trim()

    def groupingDiagnosticsConfig = config.grouping_diagnostics ?: [:]
    boolean groupingDiagnosticsRequested = groupingDiagnosticsConfig.containsKey('enabled') ? (groupingDiagnosticsConfig.enabled as boolean) : false
    if( groupingDiagnosticsRequested && !metadataPlotsEnabled ) {
        exit 1, "grouping_diagnostics.enabled requires metadata_plots.enabled to be true"
    }
    boolean groupingDiagnosticsEnabled = groupingDiagnosticsRequested
    def groupingDiagnosticsOutputDir = groupingDiagnosticsConfig.output_dir ?: 'grouping_diagnostics'
    def groupingDiagnosticsOutputDirAbs = resolveOutputRelative(groupingDiagnosticsOutputDir.toString(), outputDir)
    def groupingDiagnosticsSampleCol = groupingDiagnosticsConfig.sample_col ? groupingDiagnosticsConfig.sample_col.toString().trim() : metadataPlotsSampleCol
    def groupingDiagnosticsGroupColsRaw = groupingDiagnosticsConfig.group_cols ?: [metadataPlotsTypeCol]
    List<String> groupingDiagnosticsGroupCols = []
    if( groupingDiagnosticsGroupColsRaw instanceof List ) {
        groupingDiagnosticsGroupCols = groupingDiagnosticsGroupColsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( groupingDiagnosticsGroupColsRaw ) {
        groupingDiagnosticsGroupCols = groupingDiagnosticsGroupColsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    if( groupingDiagnosticsRequested && groupingDiagnosticsGroupCols.isEmpty() ) {
        exit 1, "grouping_diagnostics.group_cols must contain at least one metadata column when enabled"
    }
    def groupingDiagnosticsCruiseGroupColsRaw = groupingDiagnosticsConfig.cruise_level_group_cols ?: []
    List<String> groupingDiagnosticsCruiseGroupCols = groupingDiagnosticsCruiseGroupColsRaw instanceof List ?
        groupingDiagnosticsCruiseGroupColsRaw.collect { it.toString().trim() }.findAll { it } :
        groupingDiagnosticsCruiseGroupColsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    def groupingDiagnosticsCruiseCol = groupingDiagnosticsConfig.cruise_col ?: 'Cruise'
    def groupingDiagnosticsDepthCol = groupingDiagnosticsConfig.depth_col ?: 'Depth'
    def groupingDiagnosticsCruiseDepthMinPrevalence = groupingDiagnosticsConfig.cruise_depth_min_prevalence != null ? (groupingDiagnosticsConfig.cruise_depth_min_prevalence as double) : 0.5d
    def groupingDiagnosticsBaselineGroup = groupingDiagnosticsConfig.baseline_group ? groupingDiagnosticsConfig.baseline_group.toString().trim() : ''
    def groupingDiagnosticsPrimaryGroup = groupingDiagnosticsConfig.primary_group ? groupingDiagnosticsConfig.primary_group.toString().trim() : ''
    def groupingDiagnosticsAllGroupCols = (groupingDiagnosticsGroupCols + groupingDiagnosticsCruiseGroupCols).unique()
    def groupingDiagnosticsPaletteMap = extractNamedStringMap(groupingDiagnosticsConfig as Map, groupingDiagnosticsAllGroupCols, 'group_palettes', 'palette')
    def groupingDiagnosticsOrderMap = extractNamedListMap(groupingDiagnosticsConfig as Map, groupingDiagnosticsAllGroupCols, 'group_orders', 'order')
    def groupingDiagnosticsSharedPaletteConfig = config.indicspecies?.group_palettes instanceof Map ? config.indicspecies.group_palettes : [:]
    groupingDiagnosticsAllGroupCols.each { col ->
        if( !groupingDiagnosticsPaletteMap.containsKey(col) && groupingDiagnosticsSharedPaletteConfig[col] ) {
            groupingDiagnosticsPaletteMap[col] = groupingDiagnosticsSharedPaletteConfig[col]
        }
    }
    if( metadataPlotsTypeCol && metadataPlotsGroupOrder && groupingDiagnosticsGroupCols.contains(metadataPlotsTypeCol) && !groupingDiagnosticsOrderMap.containsKey(metadataPlotsTypeCol) ) {
        groupingDiagnosticsOrderMap[metadataPlotsTypeCol] = metadataPlotsGroupOrder
    }
    def groupingDiagnosticsPaletteJson = groovy.json.JsonOutput.toJson(groupingDiagnosticsPaletteMap)
    def groupingDiagnosticsOrderJson = groovy.json.JsonOutput.toJson(groupingDiagnosticsOrderMap)
    def groupingDiagnosticsMetricsRaw = groupingDiagnosticsConfig.distance_metrics ?: 'bray'
    def groupingDiagnosticsMetrics = groupingDiagnosticsMetricsRaw instanceof List ?
        groupingDiagnosticsMetricsRaw.collect { it.toString().trim().toLowerCase() }.findAll { it }.join(',') :
        groupingDiagnosticsMetricsRaw.toString().trim().toLowerCase()
    def groupingDiagnosticsTransform = groupingDiagnosticsConfig.transform ? groupingDiagnosticsConfig.transform.toString().trim().toLowerCase() : 'relative'
    def groupingDiagnosticsPermutations = groupingDiagnosticsConfig.permutations ? (groupingDiagnosticsConfig.permutations as int) : 999
    def groupingDiagnosticsRandomState = groupingDiagnosticsConfig.random_state ? (groupingDiagnosticsConfig.random_state as int) : 42
    def groupingDiagnosticsFormatsRaw = groupingDiagnosticsConfig.formats ?: 'pdf,png,svg'
    def groupingDiagnosticsFormats = groupingDiagnosticsFormatsRaw instanceof List ?
        groupingDiagnosticsFormatsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        groupingDiagnosticsFormatsRaw.toString().trim()
    def groupingDiagnosticsSoftLabelConfig = groupingDiagnosticsConfig.soft_labeling instanceof Map ? groupingDiagnosticsConfig.soft_labeling : [:]
    boolean groupingDiagnosticsSoftLabelEnabled = groupingDiagnosticsSoftLabelConfig.containsKey('enabled') ? (groupingDiagnosticsSoftLabelConfig.enabled as boolean) : false
    def groupingDiagnosticsSoftLabelK = groupingDiagnosticsSoftLabelConfig.k ? (groupingDiagnosticsSoftLabelConfig.k as int) : 7
    def groupingDiagnosticsSoftLabelTargetColsRaw = groupingDiagnosticsSoftLabelConfig.target_cols ?: (groupingDiagnosticsPrimaryGroup ? [groupingDiagnosticsPrimaryGroup] : [])
    List<String> groupingDiagnosticsSoftLabelTargetCols = groupingDiagnosticsSoftLabelTargetColsRaw instanceof List ?
        groupingDiagnosticsSoftLabelTargetColsRaw.collect { it.toString().trim() }.findAll { it } :
        groupingDiagnosticsSoftLabelTargetColsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    def groupingDiagnosticsSoftLabelExcludeRaw = groupingDiagnosticsSoftLabelConfig.exclude_labels ?: ['outlier']
    List<String> groupingDiagnosticsSoftLabelExcludeLabels = groupingDiagnosticsSoftLabelExcludeRaw instanceof List ?
        groupingDiagnosticsSoftLabelExcludeRaw.collect { it.toString().trim() }.findAll { it } :
        groupingDiagnosticsSoftLabelExcludeRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    def groupingDiagnosticsSoftLabelMinClassSamples = groupingDiagnosticsSoftLabelConfig.min_class_samples ? (groupingDiagnosticsSoftLabelConfig.min_class_samples as int) : 3
    def groupingDiagnosticsSoftLabelDistanceQuantile = groupingDiagnosticsSoftLabelConfig.distance_quantile != null ? (groupingDiagnosticsSoftLabelConfig.distance_quantile as double) : 0.95d
    boolean groupingDiagnosticsApplySoftLabels = groupingDiagnosticsSoftLabelConfig.containsKey('apply_downstream') ? (groupingDiagnosticsSoftLabelConfig.apply_downstream as boolean) : false
    def groupingDiagnosticsSoftLabelTargetCol = groupingDiagnosticsSoftLabelConfig.target_col ? groupingDiagnosticsSoftLabelConfig.target_col.toString().trim() : (groupingDiagnosticsSoftLabelTargetCols ? groupingDiagnosticsSoftLabelTargetCols[0] : '')
    def groupingDiagnosticsSoftLabelMinConfidence = groupingDiagnosticsSoftLabelConfig.min_confidence != null ? (groupingDiagnosticsSoftLabelConfig.min_confidence as double) : 0.70d
    def groupingDiagnosticsSoftLabelMinNeighborAgreement = groupingDiagnosticsSoftLabelConfig.min_neighbor_agreement != null ? (groupingDiagnosticsSoftLabelConfig.min_neighbor_agreement as double) : 0.60d
    def groupingDiagnosticsSoftLabelMinCvBalancedAccuracy = groupingDiagnosticsSoftLabelConfig.min_cv_balanced_accuracy != null ? (groupingDiagnosticsSoftLabelConfig.min_cv_balanced_accuracy as double) : 0.60d
    if( groupingDiagnosticsApplySoftLabels && (!groupingDiagnosticsEnabled || !groupingDiagnosticsSoftLabelEnabled || !groupingDiagnosticsSoftLabelTargetCol) ) {
        exit 1, "grouping_diagnostics.soft_labeling.apply_downstream requires enabled diagnostics, enabled soft labeling, and a target_col"
    }
    boolean groupingDiagnosticsRequiresMicrobialCompartments =
        (groupingDiagnosticsGroupCols + groupingDiagnosticsCruiseGroupCols).contains('microbial_compartment')
    boolean groupingDiagnosticsRunEarly = groupingDiagnosticsEnabled && !groupingDiagnosticsRequiresMicrobialCompartments
    if( groupingDiagnosticsRequiresMicrobialCompartments && !microbialCompartmentEnabled ) {
        exit 1, "grouping_diagnostics using microbial_compartment requires microbial_compartments.enabled=true"
    }
    if( groupingDiagnosticsRequiresMicrobialCompartments && groupingDiagnosticsApplySoftLabels ) {
        exit 1, "Post-clustering microbial-compartment diagnostics cannot apply inferred labels upstream"
    }
    def groupingDiagnosticsPowerConfig = groupingDiagnosticsConfig.power instanceof Map ? groupingDiagnosticsConfig.power : [:]
    boolean groupingDiagnosticsPowerEnabled = groupingDiagnosticsPowerConfig.containsKey('enabled') ? (groupingDiagnosticsPowerConfig.enabled as boolean) : false
    def groupingDiagnosticsPowerSizesRaw = groupingDiagnosticsPowerConfig.sample_sizes ?: '3,5,10,15,20'
    def groupingDiagnosticsPowerSizes = groupingDiagnosticsPowerSizesRaw instanceof List ?
        groupingDiagnosticsPowerSizesRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        groupingDiagnosticsPowerSizesRaw.toString().trim()
    def groupingDiagnosticsPowerSimulations = groupingDiagnosticsPowerConfig.simulations ? (groupingDiagnosticsPowerConfig.simulations as int) : 100
    def groupingDiagnosticsPowerPermutations = groupingDiagnosticsPowerConfig.permutations ? (groupingDiagnosticsPowerConfig.permutations as int) : 99
    def groupingDiagnosticsPowerAlpha = groupingDiagnosticsPowerConfig.alpha != null ? (groupingDiagnosticsPowerConfig.alpha as double) : 0.05d
    def groupingDiagnosticsPowerMinGroups = groupingDiagnosticsPowerConfig.min_groups ? (groupingDiagnosticsPowerConfig.min_groups as int) : 2

    def clustermapsConfig = config.clustermaps ?: [:]
    boolean clustermapsRequested = clustermapsConfig.containsKey('enabled') ? (clustermapsConfig.enabled as boolean) : false
    if( clustermapsRequested && !metadataPlotsEnabled ) {
        exit 1, "clustermaps.enabled requires metadata_plots.enabled to be true"
    }
    boolean clustermapsEnabled = clustermapsRequested
    def clustermapsOutputDir = clustermapsConfig.output_dir ?: 'clustermaps'
    def clustermapsOutputDirAbs = new File(outputDir, clustermapsOutputDir).canonicalPath
    def clustermapsMitoOutputDir = clustermapsConfig.mito_output_dir ?: 'mito/clustermaps'
    def clustermapsMitoOutputDirAbs = new File(outputDir, clustermapsMitoOutputDir).canonicalPath
    def clustermapsMitoInputPath = clustermapsConfig.mito_input ? resolveOptionalPath(clustermapsConfig.mito_input, configRoot) : new File(outputDir, 'mito/ASVs/ASV_target.mito.tsv').canonicalPath
    def clustermapsIsaFile = clustermapsConfig.isa_file ? resolveOptionalPath(clustermapsConfig.isa_file, configRoot) : null
    def clustermapsSampleCol = clustermapsConfig.sample_col ?: 'sample'
    def clustermapsSampleCodeCol = clustermapsConfig.sample_code_col ?: 'sample_code'
    def clustermapsAsvIdCol = clustermapsConfig.asv_id_col ?: 'ASV_ID'
    def clustermapsGroup1Col = clustermapsConfig.group1_col ?: 'type_group'
    def clustermapsGroup2Col = clustermapsConfig.group2_col ?: 'status'
    def clustermapsGroup3Col = clustermapsConfig.containsKey('group3_col') ? (clustermapsConfig.group3_col ?: '') : 'kit'
    def clustermapsGroup1OrderRaw = clustermapsConfig.group1_order ?: (clustermapsConfig.type_order ?: '')
    List<String> clustermapsGroup1Order = normalizePresetList(clustermapsGroup1OrderRaw, [], config.order_presets ?: [:])
    def clustermapsExcludeGroup1 = clustermapsConfig.exclude_group1 ?: (clustermapsConfig.exclude_types ?: '')
    def clustermapsGroup1Palette = clustermapsConfig.group1_palette ?: (clustermapsConfig.type_palette ?: '')
    def clustermapsGroup2Palette = clustermapsConfig.group2_palette ?: (clustermapsConfig.status_palette ?: '')
    def clustermapsGroup3Palette = clustermapsConfig.group3_palette ?: (clustermapsConfig.kit_palette ?: '')
    def clustermapsRanks = clustermapsConfig.ranks ?: 'Phylum,Class,Order,Family,Genus,Species,ASV_ID'
    def clustermapsTopN = clustermapsConfig.topN ?: 'Phylum=30,Class=30,Order=30,Family=30,Genus=30,Species=30,ASV_ID=6000'
    boolean clustermapsPlotAsvLevel = clustermapsConfig.containsKey('plot_asv_level') ? (clustermapsConfig.plot_asv_level as boolean) : true
    def clustermapsCountCol = clustermapsConfig.count_col ?: 'corr_count'
    def clustermapsIsaMinStat = clustermapsConfig.isa_min_stat != null ? (clustermapsConfig.isa_min_stat as double) : 0.6d
    def clustermapsIsaSignificanceCols = clustermapsConfig.isa_significance_cols ?: ''
    def clustermapsIsaStatCols = clustermapsConfig.isa_stat_cols ?: ''
    def clustermapsFormats = clustermapsConfig.formats ?: 'pdf,png,svg'
    def clustermapsFigWidth = clustermapsConfig.figwidth != null ? clustermapsConfig.figwidth : null
    def clustermapsRowHeight = clustermapsConfig.row_height != null ? clustermapsConfig.row_height : null
    def clustermapsMinHeight = clustermapsConfig.min_height != null ? clustermapsConfig.min_height : null
    def clustermapsMaxHeight = clustermapsConfig.max_height != null ? clustermapsConfig.max_height : null
    def clustermapsMitoSampleMode = clustermapsConfig.mito_sample_mode ?: 'auto'
    boolean clustermapsRunMito = clustermapsConfig.containsKey('run_mito') ? (clustermapsConfig.run_mito as boolean) : true
    def clustermapsIsaAutoCandidates = ([clustermapsGroup3Col, clustermapsGroup1Col, clustermapsGroup2Col] + indicspeciesGroupCols)
        .findAll { it }
        .collect { "${it}_indicator_species_summary.tsv" }
        .unique()
    def clustermapsIsaMayBeProducedInRun = clustermapsIsaFile && indicspeciesEnabled && clustermapsIsaFile == indicspeciesOutputDirAbs
    def clustermapsExternalIsaTable = null
    if( clustermapsIsaFile && !clustermapsIsaMayBeProducedInRun ) {
        File configuredIsa = new File(clustermapsIsaFile)
        if( configuredIsa.isFile() ) {
            clustermapsExternalIsaTable = configuredIsa.canonicalPath
        } else if( configuredIsa.isDirectory() ) {
            File selectedIsa = clustermapsIsaAutoCandidates
                .collect { new File(configuredIsa, it) }
                .find { it.isFile() }
            if( selectedIsa == null ) {
                selectedIsa = configuredIsa.listFiles()
                    ?.findAll { it.isFile() && it.name.endsWith('_indicator_species_summary.tsv') }
                    ?.sort { it.name }
                    ?.find()
            }
            clustermapsExternalIsaTable = selectedIsa?.canonicalPath
        }
        if( clustermapsExternalIsaTable == null ) {
            log.warn "clustermaps.isa_file did not resolve to an ISA summary table: ${clustermapsIsaFile}"
        }
    } else if( clustermapsIsaMayBeProducedInRun ) {
        log.info "CLUSTERMAPS will receive its ISA summary directly from INDICSPECIES."
    }

    def spieceasiConfig = config.spieceasi ?: [:]
    boolean spieceasiRequested = spieceasiConfig.containsKey('enabled') ? (spieceasiConfig.enabled as boolean) : false
    if( spieceasiRequested && !metadataPlotsEnabled ) {
        exit 1, "spieceasi.enabled requires metadata_plots.enabled to be true"
    }
    boolean spieceasiEnabled = spieceasiRequested
    def spieceasiOutputDir = spieceasiConfig.output_dir ?: 'spieceasi'
    def spieceasiOutputDirAbs = new File(outputDir, spieceasiOutputDir).canonicalPath
    def spieceasiPrefix = spieceasiConfig.prefix ?: 'spieceasi'
    boolean spieceasiTranspose = spieceasiConfig.containsKey('transpose') ? (spieceasiConfig.transpose as boolean) : true
    def spieceasiMinRelAbund = spieceasiConfig.min_rel_abund != null ? (spieceasiConfig.min_rel_abund as double) : 0d
    def spieceasiMinPrevalence = spieceasiConfig.min_prevalence != null ? (spieceasiConfig.min_prevalence as double) : 0.25d
    boolean spieceasiForceKeepIsaAsvs = spieceasiConfig.containsKey('force_keep_isa_asvs') ? (spieceasiConfig.force_keep_isa_asvs as boolean) : false
    if( indicspeciesRequiresMicrobialCompartments && spieceasiForceKeepIsaAsvs ) {
        exit 1, "spieceasi.force_keep_isa_asvs cannot be enabled when ISA uses microbial_compartment because that would create a circular dependency"
    }
    boolean spieceasiForceKeepAsvMagAsvs = spieceasiConfig.containsKey('force_keep_asv_mag_asvs') ? (spieceasiConfig.force_keep_asv_mag_asvs as boolean) : false
    boolean spieceasiRemoveZeroVar = spieceasiConfig.containsKey('remove_zero_var') ? (spieceasiConfig.remove_zero_var as boolean) : true
    def spieceasiMethod = spieceasiConfig.method ?: 'glasso'
    def spieceasiLambdaMinRatio = spieceasiConfig.lambda_min_ratio != null ? (spieceasiConfig.lambda_min_ratio as double) : 0.1d
    def spieceasiNlambda = spieceasiConfig.nlambda ? (spieceasiConfig.nlambda as int) : 20
    def spieceasiRepNum = spieceasiConfig.rep_num ? (spieceasiConfig.rep_num as int) : 50
    def spieceasiThresh = spieceasiConfig.thresh != null ? (spieceasiConfig.thresh as double) : 0.1d
    def spieceasiPulsarCriterion = spieceasiConfig.pulsar_criterion ? spieceasiConfig.pulsar_criterion.toString().trim().toLowerCase() : 'bstars'
    if( !['stars', 'bstars'].contains(spieceasiPulsarCriterion) ) {
        exit 1, "spieceasi.pulsar_criterion must be one of: stars, bstars"
    }
    def spieceasiNcores = spieceasiConfig.ncores ? (spieceasiConfig.ncores as int) : pipelineThreads
    def spieceasiSeed = spieceasiConfig.seed ? (spieceasiConfig.seed as int) : 10010
    def spieceasiEdgeThreshold = spieceasiConfig.edge_threshold != null ? (spieceasiConfig.edge_threshold as double) : 0.1d
    boolean spieceasiKeepNegative = spieceasiConfig.containsKey('keep_negative') ? (spieceasiConfig.keep_negative as boolean) : true
    boolean spieceasiAllPosOnly = spieceasiConfig.containsKey('all_pos_only') ? (spieceasiConfig.all_pos_only as boolean) : false
    def spieceasiLayoutIters = spieceasiConfig.layout_iters ? (spieceasiConfig.layout_iters as int) : 1000
    boolean spieceasiForceFilter = spieceasiConfig.containsKey('force_filter') ? (spieceasiConfig.force_filter as boolean) : false
    boolean spieceasiForceSpieceasi = spieceasiConfig.containsKey('force_spieceasi') ? (spieceasiConfig.force_spieceasi as boolean) : false
    boolean spieceasiForceGraphs = spieceasiConfig.containsKey('force_graphs') ? (spieceasiConfig.force_graphs as boolean) : true
    boolean networkRequested = spieceasiConfig.containsKey('network_enabled') ? (spieceasiConfig.network_enabled as boolean) : false
    if( networkRequested && !indicspeciesEnabled ) {
        exit 1, "spieceasi.network_enabled requires indicspecies.enabled to be true"
    }
    boolean networkEnabled = networkRequested && indicspeciesEnabled
    def networkGraphAllPath = spieceasiConfig.graph_pos_all ? resolveOptionalPath(spieceasiConfig.graph_pos_all, configRoot) : new File(spieceasiOutputDirAbs, "${spieceasiPrefix}_network_pos_all.graphml").canonicalPath
    def networkGraphThrPath = spieceasiConfig.graph_pos_sub ? resolveOptionalPath(spieceasiConfig.graph_pos_sub, configRoot) : new File(spieceasiOutputDirAbs, "${spieceasiPrefix}_network_pos_thr.graphml").canonicalPath
    def networkNodeFeaturesPath = spieceasiConfig.node_features ? resolveOptionalPath(spieceasiConfig.node_features, configRoot) : new File(spieceasiOutputDirAbs, "${spieceasiPrefix}_node_features.csv").canonicalPath
    if( networkEnabled && !spieceasiEnabled ) {
        [networkGraphAllPath, networkGraphThrPath, networkNodeFeaturesPath].each { p ->
            if( !new File(p).exists() ) {
                exit 1, "spieceasi.network_enabled is true while spieceasi.enabled is false, but required cached file is missing: ${p}"
            }
        }
    }
    def networkModesRaw = spieceasiConfig.network_modes
    List<String> networkModes = []
    if( networkModesRaw instanceof List ) {
        networkModes = networkModesRaw.collect { it.toString().trim() }.findAll { it }
    } else if( networkModesRaw ) {
        networkModes = networkModesRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def networkIsaOverlayGroupsRaw = spieceasiConfig.isa_overlay_groups
    List<String> networkIsaOverlayGroups = []
    if( networkIsaOverlayGroupsRaw instanceof List ) {
        networkIsaOverlayGroups = networkIsaOverlayGroupsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( networkIsaOverlayGroupsRaw ) {
        networkIsaOverlayGroups = networkIsaOverlayGroupsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    if( networkIsaOverlayGroups.isEmpty() ) {
        networkIsaOverlayGroups = indicspeciesGroupCols
    }
    networkIsaOverlayGroups = networkIsaOverlayGroups.findAll { indicspeciesGroupCols.contains(it) }
    if( networkIsaOverlayGroups.isEmpty() ) {
        networkIsaOverlayGroups = indicspeciesGroupCols
    }
    def allPosOnlyNetworkModes = [
        'degree_all',
        'abundance_all',
        'group1_isa_all',
        'group1_isa_all_labeled',
        'group1_isa_mag_all',
        'group1_isa_mag_all_labeled',
        'group2_isa_all',
        'group2_isa_all_labeled',
        'group2_isa_mag_all',
        'group2_isa_mag_all_labeled',
        'module_all',
        'module_all_labeled',
        'mag_pair_all',
        'mag_pair_all_labeled',
        'mag_pair_tax_all',
        'mag_pair_tax_all_labeled',
        'phylum_abund_all',
        'phylum_isa_all',
        'phylum_isa_all_labeled'
    ]
    if( networkIsaOverlayGroups.size() > 2 ) {
        networkIsaOverlayGroups.drop(2).eachWithIndex { groupName, offset ->
            def idx = offset + 3
            allPosOnlyNetworkModes.addAll([
                "group${idx}_isa_all",
                "group${idx}_isa_all_labeled",
                "group${idx}_isa_mag_all"
            ])
        }
    }
    if( networkModes.isEmpty() ) {
        networkModes = spieceasiAllPosOnly ? allPosOnlyNetworkModes : ['all']
    }
    if( spieceasiAllPosOnly ) {
        def modeRemap = [
            'all': null,
            'degree_sub': 'degree_all',
            'abundance_sub': 'abundance_all',
            'group1_isa': 'group1_isa_all',
            'group1_isa_labeled': 'group1_isa_all_labeled',
            'group1_isa_mag': 'group1_isa_mag_all',
            'group1_isa_mag_labeled': 'group1_isa_mag_all_labeled',
            'group1_isa_focus': 'group1_isa_focus_all',
            'group1_isa_focus_labeled': 'group1_isa_focus_all_labeled',
            'group2_isa': 'group2_isa_all',
            'group2_isa_labeled': 'group2_isa_all_labeled',
            'group2_isa_mag': 'group2_isa_mag_all',
            'group2_isa_mag_labeled': 'group2_isa_mag_all_labeled',
            'module_sub': 'module_all',
            'module_sub_labeled': 'module_all_labeled',
            'mag_pair_sub': 'mag_pair_all',
            'mag_pair_sub_labeled': 'mag_pair_all_labeled',
            'mag_pair_tax_sub': 'mag_pair_tax_all',
            'mag_pair_tax_sub_labeled': 'mag_pair_tax_all_labeled',
            'phylum_abund': 'phylum_abund_all',
            'phylum_isa': 'phylum_isa_all',
            'phylum_isa_labeled': 'phylum_isa_all_labeled'
        ]
        networkModes = networkModes.collect { mode ->
            def remapped = modeRemap.containsKey(mode) ? modeRemap[mode] : mode
            if( remapped == null ) {
                return null
            }
            def groupMode = (remapped =~ /^(group\d+_isa(?:_mag|_focus)?)(?:_labeled)?$/)
            if( groupMode.matches() ) {
                return remapped.contains('_labeled') ? "${groupMode[0][1]}_all_labeled" : "${groupMode[0][1]}_all"
            }
            return remapped
        }
            .findAll { it }
            .unique()
        if( networkModes.isEmpty() ) {
            networkModes = allPosOnlyNetworkModes
        }
    }
    def networkLayoutSeed = spieceasiConfig.layout_seed ? (spieceasiConfig.layout_seed as int) : 42
    def networkLayoutScale = spieceasiConfig.layout_scale != null ? (spieceasiConfig.layout_scale as double) : 3.0d
    def networkDegreeScale = spieceasiConfig.degree_scale != null ? (spieceasiConfig.degree_scale as double) : 80.0d
    def networkDegreeSizeMode = spieceasiConfig.degree_size_mode ? spieceasiConfig.degree_size_mode.toString().trim() : 'legacy'
    def networkDegreeMinArea = spieceasiConfig.degree_min_area != null ? (spieceasiConfig.degree_min_area as double) : 0.0d
    def networkEdgeWidthScale = spieceasiConfig.edge_width_scale != null ? (spieceasiConfig.edge_width_scale as double) : 5.0d
    def networkIsaScale = spieceasiConfig.isa_scale != null ? (spieceasiConfig.isa_scale as double) : 700.0d
    def networkAbundanceSizeMode = spieceasiConfig.abundance_size_mode ? spieceasiConfig.abundance_size_mode.toString().trim() : 'legacy'
    def networkAbundanceReference = spieceasiConfig.abundance_reference != null ? (spieceasiConfig.abundance_reference as double) : 5000.0d
    def networkAbundanceReferenceArea = spieceasiConfig.abundance_reference_area != null ? (spieceasiConfig.abundance_reference_area as double) : 80.0d
    def networkAbundanceMinArea = spieceasiConfig.abundance_min_area != null ? (spieceasiConfig.abundance_min_area as double) : 8.0d
    def networkAbundanceMaxArea = spieceasiConfig.abundance_max_area != null ? (spieceasiConfig.abundance_max_area as double) : 420.0d
    def networkAbundanceScalePower = spieceasiConfig.abundance_scale_power != null ? (spieceasiConfig.abundance_scale_power as double) : 1.6d
    def networkMaxLabels = spieceasiConfig.max_labels != null ? (spieceasiConfig.max_labels as int) : 100
    boolean networkModuleBestOnly = spieceasiConfig.containsKey('module_best_only') ? (spieceasiConfig.module_best_only as boolean) : false
    def networkModuleBestMinSize = spieceasiConfig.module_best_min_size ? (spieceasiConfig.module_best_min_size as int) : 5
    def networkModuleBestMinStability = spieceasiConfig.module_best_min_stability != null ? (spieceasiConfig.module_best_min_stability as double) : 0.7d
    def networkModuleBestTopN = spieceasiConfig.module_best_top_n != null ? (spieceasiConfig.module_best_top_n as int) : 8
    if( networkModuleBestTopN < 1 ) {
        exit 1, "spieceasi.module_best_top_n must be at least 1"
    }
    boolean networkModuleIsaOnly = spieceasiConfig.containsKey('module_isa_only') ? (spieceasiConfig.module_isa_only as boolean) : false
    boolean networkModuleSubnetworksEnabled = spieceasiConfig.containsKey('module_subnetworks_enabled') ? (spieceasiConfig.module_subnetworks_enabled as boolean) : false
    def networkModuleSubnetworkProminenceMetricsRaw = spieceasiConfig.module_subnetwork_prominence_metrics ?: ['max_relative_abundance', 'eigenvector', 'participation']
    def networkModuleSubnetworkProminenceMetrics = networkModuleSubnetworkProminenceMetricsRaw instanceof List ? networkModuleSubnetworkProminenceMetricsRaw.collect { it.toString().trim() }.findAll { it }.join(',') : networkModuleSubnetworkProminenceMetricsRaw.toString().trim()
    def networkModuleSubnetworkProminenceThreshold = spieceasiConfig.module_subnetwork_prominence_threshold != null ? (spieceasiConfig.module_subnetwork_prominence_threshold as double) : 0.75d
    def networkModuleSubnetworkLabelTopN = spieceasiConfig.module_subnetwork_label_top_n != null ? (spieceasiConfig.module_subnetwork_label_top_n as int) : 3
    def networkModuleSubnetworkProminenceMinArea = spieceasiConfig.module_subnetwork_prominence_min_area != null ? (spieceasiConfig.module_subnetwork_prominence_min_area as double) : 30.0d
    def networkModuleSubnetworkProminenceMaxArea = spieceasiConfig.module_subnetwork_prominence_max_area != null ? (spieceasiConfig.module_subnetwork_prominence_max_area as double) : 520.0d
    boolean networkModuleColorByIsa = spieceasiConfig.containsKey('module_color_by_isa') ? (spieceasiConfig.module_color_by_isa as boolean) : false
    def networkModuleIsaSource = spieceasiConfig.module_isa_source ? spieceasiConfig.module_isa_source.toString().trim() : (networkIsaOverlayGroups ? networkIsaOverlayGroups[0] : indicspeciesGroup1)
    if( networkModuleIsaSource ==~ /^group\d+$/ ) {
        def idx = networkModuleIsaSource.replaceFirst(/^group/, '') as int
        if( idx >= 1 && idx <= indicspeciesGroupCols.size() ) {
            networkModuleIsaSource = indicspeciesGroupCols[idx - 1]
        }
    }
    if( !indicspeciesGroupCols.contains(networkModuleIsaSource) ) {
        networkModuleIsaSource = indicspeciesGroupCols ? indicspeciesGroupCols[0] : 'group1'
    }
    def networkModuleIsaMinStat = spieceasiConfig.module_isa_min_stat != null ? (spieceasiConfig.module_isa_min_stat as double) : 0.25d
    def networkModuleIsaMaxQ = spieceasiConfig.module_isa_max_q != null ? (spieceasiConfig.module_isa_max_q as double) : 0.05d
    def networkMetadataPath = spieceasiConfig.metadata ? resolveOptionalPath(spieceasiConfig.metadata, configRoot) : metadataPlotsMetadataPath
    def networkColorCol = spieceasiConfig.color_col ?: indicspeciesColorCol
    def networkGroupPaletteMap = new LinkedHashMap<String,String>(indicspeciesGroupPaletteMap)
    networkGroupPaletteMap.putAll(extractNamedStringMap(spieceasiConfig as Map, indicspeciesGroupCols, 'group_palettes', 'palette'))
    def networkGroupOrderMap = new LinkedHashMap<String,List<String>>(indicspeciesGroupOrderMap)
    networkGroupOrderMap.putAll(extractNamedListMap(spieceasiConfig as Map, indicspeciesGroupCols, 'group_orders', 'order'))
    def networkFocusLabelMap = new LinkedHashMap<String,String>(indicspeciesFocusLabelMap)
    networkFocusLabelMap.putAll(extractNamedStringMap(spieceasiConfig as Map, indicspeciesGroupCols, 'focus_labels', 'focus_label'))
    def networkGroup1Palette = indicspeciesGroup1 ? (networkGroupPaletteMap[indicspeciesGroup1] ?: '') : ''
    def networkGroup2Palette = indicspeciesGroup2 ? (networkGroupPaletteMap[indicspeciesGroup2] ?: '') : ''
    List<String> networkGroup1Order = indicspeciesGroup1 ? (networkGroupOrderMap[indicspeciesGroup1] ?: []) : []
    List<String> networkGroup2Order = indicspeciesGroup2 ? (networkGroupOrderMap[indicspeciesGroup2] ?: []) : []
    def networkFocusGroup1Label = indicspeciesGroup1 ? (networkFocusLabelMap[indicspeciesGroup1] ?: '') : ''
    def networkFocusGroup2Label = indicspeciesGroup2 ? (networkFocusLabelMap[indicspeciesGroup2] ?: '') : ''
    def networkIsaOverlayGroupsCsv = networkIsaOverlayGroups.join(',')
    def networkGroupPaletteJson = groovy.json.JsonOutput.toJson(networkGroupPaletteMap)
    def networkGroupOrderJson = groovy.json.JsonOutput.toJson(networkGroupOrderMap)
    def networkFocusLabelJson = groovy.json.JsonOutput.toJson(networkFocusLabelMap)
    boolean networkModulesEnabled = networkEnabled && (spieceasiConfig.containsKey('modules_enabled') ? (spieceasiConfig.modules_enabled as boolean) : false)
    if( moduleMeasurementAssociationEnabled && !networkModulesEnabled ) {
        exit 1, "measurement_association.module_analysis_enabled requires the SPIEC-EASI network and network modules to be enabled"
    }
    def networkModuleMethodsRaw = spieceasiConfig.module_methods ?: 'leiden,louvain'
    List<String> networkModuleMethods = []
    if( networkModuleMethodsRaw instanceof List ) {
        networkModuleMethods = networkModuleMethodsRaw.collect { it.toString().trim().toLowerCase() }.findAll { it }
    } else if( networkModuleMethodsRaw ) {
        networkModuleMethods = networkModuleMethodsRaw.toString().split(/[,|]/).collect { it.trim().toLowerCase() }.findAll { it }
    }
    if( networkModuleMethods.isEmpty() ) {
        networkModuleMethods = ['leiden','louvain']
    }
    def networkModulePrimaryMethod = spieceasiConfig.module_primary_method ? spieceasiConfig.module_primary_method.toString().trim().toLowerCase() : networkModuleMethods[0]
    if( !networkModuleMethods.contains(networkModulePrimaryMethod) ) {
        networkModulePrimaryMethod = networkModuleMethods[0]
    }
    def networkModuleResolutionsRaw = spieceasiConfig.module_resolutions ?: '0.5,1.0,1.5'
    List<String> networkModuleResolutions = []
    if( networkModuleResolutionsRaw instanceof List ) {
        networkModuleResolutions = networkModuleResolutionsRaw.collect { it.toString().trim() }.findAll { it }
    } else if( networkModuleResolutionsRaw ) {
        networkModuleResolutions = networkModuleResolutionsRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    if( networkModuleResolutions.isEmpty() ) {
        networkModuleResolutions = ['1.0']
    }
    def networkModuleReps = spieceasiConfig.module_reps ? (spieceasiConfig.module_reps as int) : 25
    def networkModuleConsensusThreshold = spieceasiConfig.module_consensus_threshold != null ? (spieceasiConfig.module_consensus_threshold as double) : 0.8d
    def networkModuleSeed = spieceasiConfig.module_seed ? (spieceasiConfig.module_seed as int) : networkLayoutSeed
    def networkModulesSubPath = spieceasiConfig.modules_sub ? resolveOptionalPath(spieceasiConfig.modules_sub, configRoot) : new File(spieceasiOutputDirAbs, "${spieceasiPrefix}_modules_sub.tsv").canonicalPath
    def networkModulesAllPath = spieceasiConfig.modules_all ? resolveOptionalPath(spieceasiConfig.modules_all, configRoot) : new File(spieceasiOutputDirAbs, "${spieceasiPrefix}_modules_all.tsv").canonicalPath

    def genomeCooccurrenceConfig = config.genome_cooccurrence ?: [:]
    boolean genomeCooccurrenceEnabled = genomeCooccurrenceConfig.containsKey('enabled') ? (genomeCooccurrenceConfig.enabled as boolean) : false
    def genomeCooccurrenceOutputDir = genomeCooccurrenceConfig.output_dir ?: 'genome_cooccurrence'
    def genomeCooccurrenceOutputDirAbs = resolveOutputRelative(genomeCooccurrenceOutputDir.toString(), outputDir)
    def genomeCooccurrenceGenomeQc = genomeCooccurrenceConfig.genome_qc ? resolveOptionalPath(genomeCooccurrenceConfig.genome_qc.toString(), configRoot) : null
    def genomeCooccurrenceMetagenome = genomeCooccurrenceConfig.metagenome_abundance ? resolveOptionalPath(genomeCooccurrenceConfig.metagenome_abundance.toString(), configRoot) : null
    def genomeCooccurrenceMetatranscriptome = genomeCooccurrenceConfig.metatranscriptome_abundance ? resolveOptionalPath(genomeCooccurrenceConfig.metatranscriptome_abundance.toString(), configRoot) : null
    def genomeCooccurrenceMetadata = genomeCooccurrenceConfig.metadata ? resolveOptionalPath(genomeCooccurrenceConfig.metadata.toString(), configRoot) : null
    def genomeCooccurrenceSampleCol = genomeCooccurrenceConfig.sample_col ?: 'sample_id'
    def genomeCooccurrenceGenomeCol = genomeCooccurrenceConfig.genome_col ?: 'genome_id'
    def genomeCooccurrenceMetagenomeValueCol = genomeCooccurrenceConfig.metagenome_value_col ?: 'read_count'
    def genomeCooccurrenceMetatranscriptomeValueCol = genomeCooccurrenceConfig.metatranscriptome_value_col ?: 'read_count'
    def genomeCooccurrenceInputScale = genomeCooccurrenceConfig.input_scale ?: 'raw_counts'
    def genomeCooccurrenceInputFeatureLevel = genomeCooccurrenceConfig.input_feature_level ?: 'species'
    def genomeCooccurrenceExpectedSpecies = genomeCooccurrenceConfig.expected_species != null ? (genomeCooccurrenceConfig.expected_species as int) : 24
    def genomeCooccurrenceClosureTotal = genomeCooccurrenceConfig.closure_total != null ? (genomeCooccurrenceConfig.closure_total as double) : 1000000.0d
    def genomeCooccurrenceMinimumReadCount = genomeCooccurrenceConfig.minimum_read_count != null ? (genomeCooccurrenceConfig.minimum_read_count as double) : 0.0d
    def genomeCooccurrenceMetadataSampleCol = genomeCooccurrenceConfig.metadata_sample_col ?: 'sample'
    def genomeCooccurrenceCruiseCol = genomeCooccurrenceConfig.cruise_col ?: 'Cruise'
    def genomeCooccurrenceSeasonCol = genomeCooccurrenceConfig.season_col ?: 'Season'
    def genomeCooccurrenceDepthCol = genomeCooccurrenceConfig.depth_col ?: 'Depth'
    def genomeCooccurrenceMonthCol = genomeCooccurrenceConfig.month_col ?: 'Month'
    def genomeCooccurrenceSampleIdRegex = genomeCooccurrenceConfig.sample_id_regex ?: '^SI(?P<cruise>[0-9]+)_(?P<depth>[0-9]+(?:[.][0-9]+)?)m$'
    def genomeCooccurrenceMinRelAbund = genomeCooccurrenceConfig.min_rel_abund != null ? (genomeCooccurrenceConfig.min_rel_abund as double) : 0.0d
    def genomeCooccurrenceMinPrevalence = genomeCooccurrenceConfig.min_prevalence != null ? (genomeCooccurrenceConfig.min_prevalence as double) : 0.05d
    def genomeCooccurrenceZeroReplacementFraction = genomeCooccurrenceConfig.zero_replacement_fraction != null ? (genomeCooccurrenceConfig.zero_replacement_fraction as double) : 0.65d
    def genomeCooccurrenceMinAbsRho = genomeCooccurrenceConfig.min_abs_rho != null ? (genomeCooccurrenceConfig.min_abs_rho as double) : 0.30d
    def genomeCooccurrenceBootstrapIterations = genomeCooccurrenceConfig.bootstrap_iterations != null ? (genomeCooccurrenceConfig.bootstrap_iterations as int) : 1000
    def genomeCooccurrenceMinBootstrapRecovery = genomeCooccurrenceConfig.min_bootstrap_recovery != null ? (genomeCooccurrenceConfig.min_bootstrap_recovery as double) : 0.80d
    def genomeCooccurrenceMinSignConsistency = genomeCooccurrenceConfig.min_sign_consistency != null ? (genomeCooccurrenceConfig.min_sign_consistency as double) : 0.90d
    def genomeCooccurrencePermutations = genomeCooccurrenceConfig.permutations != null ? (genomeCooccurrenceConfig.permutations as int) : 1000
    def genomeCooccurrencePermutationStrata = genomeCooccurrenceConfig.permutation_strata instanceof List ? genomeCooccurrenceConfig.permutation_strata.join(',') : (genomeCooccurrenceConfig.permutation_strata ?: 'Season,Depth')
    def genomeCooccurrenceMaxQ = genomeCooccurrenceConfig.max_q != null ? (genomeCooccurrenceConfig.max_q as double) : 0.05d
    def genomeCooccurrenceSeed = genomeCooccurrenceConfig.seed != null ? (genomeCooccurrenceConfig.seed as int) : 10010
    if( genomeCooccurrenceEnabled ) {
        [
            'genome_cooccurrence.genome_qc': genomeCooccurrenceGenomeQc,
            'genome_cooccurrence.metagenome_abundance': genomeCooccurrenceMetagenome,
            'genome_cooccurrence.metatranscriptome_abundance': genomeCooccurrenceMetatranscriptome,
            'genome_cooccurrence.metadata': genomeCooccurrenceMetadata,
        ].each { label, pathValue ->
            if( !pathValue || !new File(pathValue).isFile() ) {
                exit 1, "${label} file not found: ${pathValue}"
            }
        }
    }

    def masterSummaryConfig = config.master_summary ?: [:]
    boolean masterSummaryEnabled = masterSummaryConfig.containsKey('enabled') ? (masterSummaryConfig.enabled as boolean) : false
    if( masterSummaryEnabled && !metadataPlotsEnabled ) {
        exit 1, "master_summary.enabled requires metadata_plots.enabled to be true"
    }
    def masterSummaryOutputDir = masterSummaryConfig.output_dir ?: 'summary/tables'
    def masterSummaryOutputDirAbs = resolveOutputRelative(masterSummaryOutputDir.toString(), outputDir)
    def masterSummaryClustermapsDir = masterSummaryConfig.clustermaps_dir ?: 'clustermaps'
    def masterSummaryClustermapsDirAbs = resolveOutputRelative(masterSummaryClustermapsDir.toString(), outputDir)
    def masterSummaryIndicspeciesDir = masterSummaryConfig.indicspecies_dir ?: 'indicspecies'
    def masterSummaryIndicspeciesDirAbs = resolveOutputRelative(masterSummaryIndicspeciesDir.toString(), outputDir)
    def masterSummarySpieceasiDir = masterSummaryConfig.spieceasi_dir ?: 'spieceasi'
    def masterSummarySpieceasiDirAbs = resolveOutputRelative(masterSummarySpieceasiDir.toString(), outputDir)
    def masterSummaryAsvMagDir = masterSummaryConfig.asv_mag_dir ?: 'asv_mag_link'
    def masterSummaryAsvMagDirAbs = resolveOutputRelative(masterSummaryAsvMagDir.toString(), outputDir)
    def masterSummaryWhitelistRaw = masterSummaryConfig.whitelist
    List<String> masterSummaryWhitelist = []
    if( masterSummaryWhitelistRaw instanceof List ) {
        masterSummaryWhitelist = masterSummaryWhitelistRaw.collect { it.toString().trim() }.findAll { it }
    } else if( masterSummaryWhitelistRaw ) {
        masterSummaryWhitelist = masterSummaryWhitelistRaw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
    }
    def masterSummaryWhitelistCsv = masterSummaryWhitelist ? masterSummaryWhitelist.join(',') : ''
    def masterSummaryMaxDirectCols = masterSummaryConfig.max_direct_cols ? (masterSummaryConfig.max_direct_cols as int) : 300

    def asvMagLinkConfig = config.asv_mag_link ?: [:]
    boolean asvMagLinkEnabled = asvMagLinkConfig.containsKey('enabled') ? (asvMagLinkConfig.enabled as boolean) : false
    def asvMagLinkMasterTsv = asvMagLinkConfig.master_tsv ? resolveOptionalPath(asvMagLinkConfig.master_tsv, configRoot) : null
    def asvMagLinkBarrnapDir = asvMagLinkConfig.barrnap_dir ? resolveOptionalPath(asvMagLinkConfig.barrnap_dir, configRoot) : null
    def asvMagLinkGenomeDir = asvMagLinkConfig.genome_fasta_dir ? resolveOptionalPath(asvMagLinkConfig.genome_fasta_dir, configRoot) : null
    def asvMagLinkGenomeQcDir = asvMagLinkConfig.genome_qc_dir ? resolveOptionalPath(asvMagLinkConfig.genome_qc_dir, configRoot) : null
    boolean asvMagLinkAutoBarrnap = asvMagLinkConfig.containsKey('auto_barrnap') ? (asvMagLinkConfig.auto_barrnap as boolean) : false
    def asvMagLinkGenomeQcDirsRaw = asvMagLinkConfig.genome_qc_dirs
    def asvMagLinkIdTokenIndexesRaw = asvMagLinkConfig.id_token_indexes
    List asvMagLinkGenomeQcDirs = []
    List asvMagLinkIdTokenIndexes = []
    if( asvMagLinkGenomeQcDirsRaw instanceof List ) {
        asvMagLinkGenomeQcDirs = asvMagLinkGenomeQcDirsRaw.collect { resolveOptionalPath(it, configRoot) }.findAll { it }
    } else if( asvMagLinkGenomeQcDirsRaw ) {
        asvMagLinkGenomeQcDirs = asvMagLinkGenomeQcDirsRaw.toString().split(/[,|]/).collect { resolveOptionalPath(it.trim(), configRoot) }.findAll { it }
    }
    if( asvMagLinkIdTokenIndexesRaw instanceof List ) {
        asvMagLinkIdTokenIndexes = asvMagLinkIdTokenIndexesRaw.collect { it == null || it.toString().trim() == '' ? -1 : (it as int) }
    } else if( asvMagLinkIdTokenIndexesRaw ) {
        asvMagLinkIdTokenIndexes = asvMagLinkIdTokenIndexesRaw.toString().split(/[,|]/).collect { token ->
            def trimmed = token.trim()
            trimmed ? (trimmed as int) : -1
        }
    }
    def asvMagLinkOutputDir = asvMagLinkConfig.output_dir ?: 'asv_mag_link'
    def asvMagLinkOutputDirAbs = resolveOutputRelative(asvMagLinkOutputDir.toString(), outputDir)
    def asvMagLinkThreads = asvMagLinkConfig.threads ? (asvMagLinkConfig.threads as int) : pipelineThreads
    def asvMagLinkMinPident = asvMagLinkConfig.min_pident != null ? (asvMagLinkConfig.min_pident as double) : 97.0d
    def asvMagLinkMinQcov = asvMagLinkConfig.min_qcov != null ? (asvMagLinkConfig.min_qcov as double) : 90.0d
    def asvMagLinkAsvTaxonomyMinConfidence = asvMagLinkConfig.asv_taxonomy_min_confidence != null ?
        (asvMagLinkConfig.asv_taxonomy_min_confidence as double) : null
    if( asvMagLinkAsvTaxonomyMinConfidence != null &&
        (asvMagLinkAsvTaxonomyMinConfidence < 0d || asvMagLinkAsvTaxonomyMinConfidence > 1d) ) {
        exit 1, "asv_mag_link.asv_taxonomy_min_confidence must be between 0 and 1"
    }
    def asvMagLinkMinCompleteness = asvMagLinkConfig.min_completeness != null ? (asvMagLinkConfig.min_completeness as double) : null
    def asvMagLinkMaxContamination = asvMagLinkConfig.max_contamination != null ? (asvMagLinkConfig.max_contamination as double) : null
    def asvMagLinkGuncAssessmentValue = asvMagLinkConfig.gunc_assessment_value ?
        asvMagLinkConfig.gunc_assessment_value.toString().trim() : ''
    boolean asvMagLinkRequireSpeciesAssignment = asvMagLinkConfig.containsKey('require_species_assignment') ?
        (asvMagLinkConfig.require_species_assignment as boolean) : false
    def asvMagLinkMinRrnaMarkerCount = asvMagLinkConfig.min_rrna_marker_count != null ?
        (asvMagLinkConfig.min_rrna_marker_count as double) : null
    def asvMagLinkTopN = asvMagLinkConfig.top_n ? (asvMagLinkConfig.top_n as int) : 5
    def asvMagLinkPlotTopN = asvMagLinkConfig.plot_top_n ? (asvMagLinkConfig.plot_top_n as int) : 20
    if( asvMagLinkEnabled && asvMagLinkMasterTsv && (asvMagLinkBarrnapDir || asvMagLinkGenomeDir || asvMagLinkGenomeQcDir || asvMagLinkGenomeQcDirs) ) {
        exit 1, "asv_mag_link.master_tsv cannot be combined with genome_qc_dir, genome_qc_dirs, barrnap_dir, or genome_fasta_dir"
    }
    if( asvMagLinkEnabled && !asvMagLinkMasterTsv && !asvMagLinkBarrnapDir && !asvMagLinkGenomeQcDir && !asvMagLinkGenomeQcDirs ) {
        exit 1, "asv_mag_link.enabled requires asv_mag_link.master_tsv, genome_qc_dir, genome_qc_dirs, or barrnap_dir"
    }

    def asvMagNetworkConfig = config.asv_mag_network ?: [:]
    boolean asvMagNetworkRequested = asvMagNetworkConfig.containsKey('enabled') ? (asvMagNetworkConfig.enabled as boolean) : false
    if( asvMagNetworkRequested && !networkEnabled ) {
        exit 1, "asv_mag_network.enabled requires network.enabled"
    }
    if( asvMagNetworkRequested && !asvMagLinkEnabled ) {
        exit 1, "asv_mag_network.enabled requires asv_mag_link.enabled"
    }
    boolean asvMagNetworkEnabled = asvMagNetworkRequested
    if( genomeCooccurrenceEnabled && (!asvMagNetworkEnabled || !networkModulesEnabled || !titanEnabled || !indicspeciesEnabled) ) {
        exit 1, "genome_cooccurrence requires ASV-MAG network mappings, ASV ecological modules, TITAN, and indicator-species analyses"
    }
    def asvMagNetworkOutputDir = asvMagNetworkConfig.output_dir ?: 'asv_mag_network'
    def asvMagNetworkOutputDirAbs = resolveOutputRelative(asvMagNetworkOutputDir.toString(), outputDir)
    def asvMagNetworkPrefix = asvMagNetworkConfig.prefix ?: 'asv_mag_network'
    def asvMagNetworkAnchorTopN = asvMagNetworkConfig.anchor_top_n != null ? (asvMagNetworkConfig.anchor_top_n as int) : 1
    if( asvMagNetworkAnchorTopN < 1 ) {
        exit 1, "asv_mag_network.anchor_top_n must be at least 1"
    }
    def asvMagNetworkGraphVariant = asvMagNetworkConfig.graph_variant ? asvMagNetworkConfig.graph_variant.toString().trim().toLowerCase() : (spieceasiAllPosOnly ? 'all' : 'thresholded')
    if( !['all', 'thresholded'].contains(asvMagNetworkGraphVariant) ) {
        exit 1, "asv_mag_network.graph_variant must be one of: all, thresholded"
    }
    def asvMagNetworkMinPident = asvMagNetworkConfig.min_pident != null ? (asvMagNetworkConfig.min_pident as double) : 99.5d
    def asvMagNetworkMinQcov = asvMagNetworkConfig.min_qcov != null ? (asvMagNetworkConfig.min_qcov as double) : 100.0d
    def asvMagNetworkAsvTaxonomySource = asvMagNetworkConfig.asv_taxonomy_source ? asvMagNetworkConfig.asv_taxonomy_source.toString().trim() : 'ncbi'
    def asvMagNetworkMagTaxonomySource = asvMagNetworkConfig.mag_taxonomy_source ? asvMagNetworkConfig.mag_taxonomy_source.toString().trim() : 'gtdb'
    def asvMagNetworkMagAbundance = asvMagNetworkConfig.mag_abundance ? resolveOptionalPath(asvMagNetworkConfig.mag_abundance, configRoot) : null
    asvMagNetworkCruiseMetadata = asvMagNetworkConfig.cruise_metadata ? resolveOptionalPath(asvMagNetworkConfig.cruise_metadata, configRoot) : null
    def asvMagNetworkMagIdMode = asvMagNetworkConfig.mag_id_mode ? asvMagNetworkConfig.mag_id_mode.toString().trim().toLowerCase() : 'exact'
    if( !['exact', 'suffix_after_double_underscore'].contains(asvMagNetworkMagIdMode) ) {
        exit 1, "asv_mag_network.mag_id_mode must be one of: exact, suffix_after_double_underscore"
    }
    def asvMagNetworkMagAbundanceFormat = asvMagNetworkConfig.mag_abundance_format ? asvMagNetworkConfig.mag_abundance_format.toString().trim().toLowerCase() : 'auto'
    if( !['auto', 'long', 'wide'].contains(asvMagNetworkMagAbundanceFormat) ) {
        exit 1, "asv_mag_network.mag_abundance_format must be one of: auto, long, wide"
    }
    def asvMagNetworkMagAbundanceGenomeCol = asvMagNetworkConfig.mag_abundance_genome_col ? asvMagNetworkConfig.mag_abundance_genome_col.toString().trim() : 'genome_id'
    def asvMagNetworkMagAbundanceSampleCol = asvMagNetworkConfig.mag_abundance_sample_col ? asvMagNetworkConfig.mag_abundance_sample_col.toString().trim() : 'sample_id'
    def asvMagNetworkMagAbundanceValueCol = asvMagNetworkConfig.mag_abundance_value_col ? asvMagNetworkConfig.mag_abundance_value_col.toString().trim() : 'read_count'
    def asvMagNetworkMagAbundanceSeqkit = asvMagNetworkConfig.mag_abundance_seqkit ? resolveOptionalPath(asvMagNetworkConfig.mag_abundance_seqkit, configRoot) : null
    def asvMagNetworkMagAbundanceSeqkitFileCol = asvMagNetworkConfig.mag_abundance_seqkit_file_col ? asvMagNetworkConfig.mag_abundance_seqkit_file_col.toString().trim() : 'file'
    def asvMagNetworkMagAbundanceSeqkitCountCol = asvMagNetworkConfig.mag_abundance_seqkit_count_col ? asvMagNetworkConfig.mag_abundance_seqkit_count_col.toString().trim() : 'num_seqs'
    def asvMagNetworkMagTranscriptAbundance = asvMagNetworkConfig.mag_transcript_abundance ? resolveOptionalPath(asvMagNetworkConfig.mag_transcript_abundance, configRoot) : null
    def asvMagNetworkMagTranscriptAbundanceFormat = asvMagNetworkConfig.mag_transcript_abundance_format ? asvMagNetworkConfig.mag_transcript_abundance_format.toString().trim().toLowerCase() : 'auto'
    if( !['auto', 'long', 'wide'].contains(asvMagNetworkMagTranscriptAbundanceFormat) ) {
        exit 1, "asv_mag_network.mag_transcript_abundance_format must be one of: auto, long, wide"
    }
    def asvMagNetworkMagTranscriptAbundanceGenomeCol = asvMagNetworkConfig.mag_transcript_abundance_genome_col ? asvMagNetworkConfig.mag_transcript_abundance_genome_col.toString().trim() : 'genome_id'
    def asvMagNetworkMagTranscriptAbundanceSampleCol = asvMagNetworkConfig.mag_transcript_abundance_sample_col ? asvMagNetworkConfig.mag_transcript_abundance_sample_col.toString().trim() : 'sample_id'
    def asvMagNetworkMagTranscriptAbundanceValueCol = asvMagNetworkConfig.mag_transcript_abundance_value_col ? asvMagNetworkConfig.mag_transcript_abundance_value_col.toString().trim() : 'read_count'
    def asvMagNetworkMagTranscriptAbundanceSeqkit = asvMagNetworkConfig.mag_transcript_abundance_seqkit ? resolveOptionalPath(asvMagNetworkConfig.mag_transcript_abundance_seqkit, configRoot) : null
    def asvMagNetworkMagTranscriptAbundanceSeqkitFileCol = asvMagNetworkConfig.mag_transcript_abundance_seqkit_file_col ? asvMagNetworkConfig.mag_transcript_abundance_seqkit_file_col.toString().trim() : 'file'
    def asvMagNetworkMagTranscriptAbundanceSeqkitCountCol = asvMagNetworkConfig.mag_transcript_abundance_seqkit_count_col ? asvMagNetworkConfig.mag_transcript_abundance_seqkit_count_col.toString().trim() : 'num_seqs'
    def asvMagNetworkMagTranscriptAbundanceNormalization = asvMagNetworkConfig.mag_transcript_abundance_normalization ? asvMagNetworkConfig.mag_transcript_abundance_normalization.toString().trim().toLowerCase() : 'auto'
    if( !['auto', 'input_fragment_fpm', 'provided_fpkm', 'provided_tpm', 'median_ratio', 'raw', 'relative'].contains(asvMagNetworkMagTranscriptAbundanceNormalization) ) {
        exit 1, "asv_mag_network.mag_transcript_abundance_normalization must be one of: auto, input_fragment_fpm, provided_fpkm, provided_tpm, median_ratio, raw, relative"
    }
    def asvMagNetworkMinSharedSamples = asvMagNetworkConfig.min_shared_samples ? (asvMagNetworkConfig.min_shared_samples as int) : 5
    def asvMagNetworkMagKnn = asvMagNetworkConfig.mag_knn ? (asvMagNetworkConfig.mag_knn as int) : 3
    if( asvMagNetworkMagKnn < 1 ) {
        exit 1, "asv_mag_network.mag_knn must be at least 1"
    }
    def asvMagNetworkAbundanceTransform = asvMagNetworkConfig.abundance_transform ? asvMagNetworkConfig.abundance_transform.toString().trim().toLowerCase() : 'log1p'
    if( !['none', 'log1p'].contains(asvMagNetworkAbundanceTransform) ) {
        exit 1, "asv_mag_network.abundance_transform must be one of: none, log1p"
    }
    def asvMagNetworkMagAbundanceNormalization = asvMagNetworkConfig.mag_abundance_normalization ? asvMagNetworkConfig.mag_abundance_normalization.toString().trim().toLowerCase() : 'auto'
    if( !['auto', 'input_fragment_fpm', 'provided_fpkm', 'provided_tpm', 'median_ratio', 'raw', 'relative'].contains(asvMagNetworkMagAbundanceNormalization) ) {
        exit 1, "asv_mag_network.mag_abundance_normalization must be one of: auto, input_fragment_fpm, provided_fpkm, provided_tpm, median_ratio, raw, relative"
    }
    if( asvMagNetworkMagAbundanceSeqkit && !asvMagNetworkMagAbundance ) {
        exit 1, "asv_mag_network.mag_abundance_seqkit requires mag_abundance"
    }
    if( asvMagNetworkMagTranscriptAbundanceSeqkit && !asvMagNetworkMagTranscriptAbundance ) {
        exit 1, "asv_mag_network.mag_transcript_abundance_seqkit requires mag_transcript_abundance"
    }
    if( asvMagNetworkMagAbundance && (asvMagNetworkMagAbundanceNormalization in ['auto', 'input_fragment_fpm']) && !asvMagNetworkMagAbundanceSeqkit ) {
        exit 1, "asv_mag_network.mag_abundance_normalization=${asvMagNetworkMagAbundanceNormalization} requires mag_abundance_seqkit; select median_ratio explicitly to use the legacy normalization"
    }
    if( asvMagNetworkMagTranscriptAbundance && (asvMagNetworkMagTranscriptAbundanceNormalization in ['auto', 'input_fragment_fpm']) && !asvMagNetworkMagTranscriptAbundanceSeqkit ) {
        exit 1, "asv_mag_network.mag_transcript_abundance_normalization=${asvMagNetworkMagTranscriptAbundanceNormalization} requires mag_transcript_abundance_seqkit; select median_ratio explicitly to use the legacy normalization"
    }
    [
        'mag_abundance': asvMagNetworkMagAbundance,
        'mag_abundance_seqkit': asvMagNetworkMagAbundanceSeqkit,
        'mag_transcript_abundance': asvMagNetworkMagTranscriptAbundance,
        'mag_transcript_abundance_seqkit': asvMagNetworkMagTranscriptAbundanceSeqkit,
    ].each { label, candidate ->
        if( candidate && !new File(candidate).isFile() ) {
            exit 1, "asv_mag_network.${label} not found: ${candidate}"
        }
    }
    def asvMagNetworkFunctionalModuleMinFraction = asvMagNetworkConfig.functional_module_min_fraction != null ? (asvMagNetworkConfig.functional_module_min_fraction as double) : 0.5d
    def asvMagNetworkAmbiguityTargetModulesRaw = asvMagNetworkConfig.ambiguity_target_modules ?: ['auto_n_s']
    def asvMagNetworkAmbiguityTargetModules = asvMagNetworkAmbiguityTargetModulesRaw instanceof List ? asvMagNetworkAmbiguityTargetModulesRaw.join(',') : asvMagNetworkAmbiguityTargetModulesRaw.toString()
    def asvMagNetworkAmbiguityFormats = asvMagNetworkConfig.ambiguity_formats ? asvMagNetworkConfig.ambiguity_formats.toString().trim() : 'pdf,png,svg'
    def asvMagNetworkFontFamily = asvMagNetworkConfig.font_family ? asvMagNetworkConfig.font_family.toString().trim() : 'Times New Roman'
    def asvMagNetworkBiochemicalGroupingsRaw = asvMagNetworkConfig.biochemical_groupings ?: ['o2_subcompartment_final', 'cruise_group']
    def asvMagNetworkBiochemicalGroupings = asvMagNetworkBiochemicalGroupingsRaw instanceof List ?
        asvMagNetworkBiochemicalGroupingsRaw.join(',') :
        asvMagNetworkBiochemicalGroupingsRaw.toString()
    def asvMagNetworkIsaQThreshold = asvMagNetworkConfig.isa_q_threshold != null ?
        (asvMagNetworkConfig.isa_q_threshold as double) : 0.05d
    def asvMagNetworkGroupingSampleCol = asvMagNetworkConfig.grouping_diagnostic_sample_col ?: 'sample'
    def asvMagNetworkGroupingCruiseCol = asvMagNetworkConfig.grouping_diagnostic_cruise_col ?: 'Cruise'
    def asvMagNetworkGroupingGroupsRaw = asvMagNetworkConfig.grouping_diagnostic_groups ?: ['o2_subcompartment_final', 'cruise_group']
    def asvMagNetworkGroupingGroups = asvMagNetworkGroupingGroupsRaw instanceof List ?
        asvMagNetworkGroupingGroupsRaw.join(',') :
        asvMagNetworkGroupingGroupsRaw.toString()
    def asvMagNetworkGroupingCruiseGroupsRaw = asvMagNetworkConfig.grouping_diagnostic_cruise_level_groups ?: ['Season', 'cruise_group']
    def asvMagNetworkGroupingCruiseGroups = asvMagNetworkGroupingCruiseGroupsRaw instanceof List ?
        asvMagNetworkGroupingCruiseGroupsRaw.join(',') :
        asvMagNetworkGroupingCruiseGroupsRaw.toString()
    def asvMagNetworkGroupingPermutations = asvMagNetworkConfig.grouping_diagnostic_permutations != null ?
        (asvMagNetworkConfig.grouping_diagnostic_permutations as int) : 999
    def asvMagNetworkFunctionalRaw = asvMagNetworkConfig.functional_annotations ?: []
    List asvMagNetworkFunctionalAnnotations = []
    if( asvMagNetworkFunctionalRaw instanceof List ) {
        asvMagNetworkFunctionalAnnotations = asvMagNetworkFunctionalRaw.collect { resolveOptionalPath(it, configRoot) }.findAll { it }
    } else if( asvMagNetworkFunctionalRaw ) {
        asvMagNetworkFunctionalAnnotations = asvMagNetworkFunctionalRaw.toString().split(/[,|]/).collect { resolveOptionalPath(it.trim(), configRoot) }.findAll { it }
    }

    def asvMagCurtainsConfig = config.asv_mag_curtains ?: [:]
    boolean asvMagCurtainsEnabled = asvMagCurtainsConfig.containsKey('enabled') ? (asvMagCurtainsConfig.enabled as boolean) : false
    if( asvMagCurtainsEnabled && (!asvMagNetworkEnabled || !microbialStateInterpretationEnabled) ) {
        exit 1, "asv_mag_curtains.enabled requires asv_mag_network.enabled and microbial_state_interpretation.enabled"
    }
    if( asvMagCurtainsEnabled && (!asvMagNetworkMagAbundance || !asvMagNetworkMagTranscriptAbundance) ) {
        exit 1, "asv_mag_curtains.enabled requires both metagenome and metatranscriptome abundance inputs"
    }
    def asvMagCurtainsOutputDir = asvMagCurtainsConfig.output_dir ?: 'asv_mag_curtains'
    def asvMagCurtainsOutputDirAbs = resolveOutputRelative(asvMagCurtainsOutputDir.toString(), outputDir)
    def asvMagCurtainsFormatsRaw = asvMagCurtainsConfig.formats ?: 'pdf,png,svg'
    def asvMagCurtainsFormats = asvMagCurtainsFormatsRaw instanceof List ? asvMagCurtainsFormatsRaw.join(',') : asvMagCurtainsFormatsRaw.toString()
    double asvMagCurtainsMaximumDepth = asvMagCurtainsConfig.maximum_depth != null ? (asvMagCurtainsConfig.maximum_depth as double) : microbialStateInterpretationMaximumDepth
    double asvMagCurtainsMinimumPointArea = asvMagCurtainsConfig.minimum_point_area != null ? (asvMagCurtainsConfig.minimum_point_area as double) : 1.0d
    double asvMagCurtainsMaximumPointArea = asvMagCurtainsConfig.maximum_point_area != null ? (asvMagCurtainsConfig.maximum_point_area as double) : 320.0d
    double asvMagCurtainsAsvSamplePointArea = asvMagCurtainsConfig.asv_sample_point_area != null ? (asvMagCurtainsConfig.asv_sample_point_area as double) : 6.0d
    double asvMagCurtainsContrastPseudocountFpm = asvMagCurtainsConfig.contrast_pseudocount_fpm != null ? (asvMagCurtainsConfig.contrast_pseudocount_fpm as double) : 1.0d
    double asvMagCurtainsContrastFoldThreshold = asvMagCurtainsConfig.contrast_fold_threshold != null ? (asvMagCurtainsConfig.contrast_fold_threshold as double) : 2.0d
    def asvMagCurtainsRenewalEvents = asvMagCurtainsConfig.renewal_events ?
        resolveOptionalPath(asvMagCurtainsConfig.renewal_events.toString(), configRoot) : null
    def asvMagCurtainsRenewalDateCol = (asvMagCurtainsConfig.renewal_date_col ?: 'start_date').toString()
    if( asvMagCurtainsEnabled && (asvMagCurtainsMaximumDepth <= 0 || asvMagCurtainsMinimumPointArea <= 0 ||
        asvMagCurtainsMaximumPointArea <= asvMagCurtainsMinimumPointArea || asvMagCurtainsAsvSamplePointArea <= 0 ||
        asvMagCurtainsContrastPseudocountFpm <= 0 || asvMagCurtainsContrastFoldThreshold <= 1) ) {
        exit 1, "asv_mag_curtains depth and point-area settings are invalid"
    }

    def groupGuildFunctionConfig = config.group_guild_function ?: [:]
    boolean groupGuildFunctionEnabled = groupGuildFunctionConfig.containsKey('enabled') ? (groupGuildFunctionConfig.enabled as boolean) : false
    if( groupGuildFunctionEnabled && (!asvMagNetworkEnabled || !networkModulesEnabled || !indicspeciesEnabled) ) {
        exit 1, "group_guild_function.enabled requires asv_mag_network, network_modules, and indicspecies"
    }
    def groupGuildFunctionOutputDirAbs = resolveOutputRelative((groupGuildFunctionConfig.output_dir ?: 'group_guild_function').toString(), outputDir)
    def groupGuildFunctionGroupingsRaw = groupGuildFunctionConfig.groupings ?: ['o2_subcompartment_final','cruise_group','renewal_phase','o2_compartment','gmm_component']
    def groupGuildFunctionGroupings = groupGuildFunctionGroupingsRaw instanceof List ? groupGuildFunctionGroupingsRaw.join(',') : groupGuildFunctionGroupingsRaw.toString()
    def groupGuildFunctionPrimaryRaw = groupGuildFunctionConfig.primary_groupings ?: ['o2_compartment','gmm_component','o2_subcompartment_final']
    def groupGuildFunctionPrimary = groupGuildFunctionPrimaryRaw instanceof List ? groupGuildFunctionPrimaryRaw.join(',') : groupGuildFunctionPrimaryRaw.toString()
    def groupGuildFunctionTargetsRaw = groupGuildFunctionConfig.target_modules ?: ['auto_n_s']
    def groupGuildFunctionTargets = groupGuildFunctionTargetsRaw instanceof List ? groupGuildFunctionTargetsRaw.join(',') : groupGuildFunctionTargetsRaw.toString()
    def groupGuildFunctionSampleCol = groupGuildFunctionConfig.sample_col ?: metadataPlotsSampleCol
    def groupGuildFunctionCruiseCol = groupGuildFunctionConfig.cruise_col ?: 'Cruise'
    def groupGuildFunctionPermutations = groupGuildFunctionConfig.permutations != null ? (groupGuildFunctionConfig.permutations as int) : 999
    def groupGuildFunctionSeed = groupGuildFunctionConfig.seed != null ? (groupGuildFunctionConfig.seed as int) : 42
    def groupGuildFunctionFormats = groupGuildFunctionConfig.formats ?: 'pdf,png,svg'
    def groupGuildFunctionFontFamily = groupGuildFunctionConfig.font_family ?: 'Times New Roman'
    boolean groupGuildFunctionStrictFont = groupGuildFunctionConfig.containsKey('strict_font') ? (groupGuildFunctionConfig.strict_font as boolean) : false
    def groupGuildFunctionMaxModules = groupGuildFunctionConfig.max_modules != null ? (groupGuildFunctionConfig.max_modules as int) : 8
    def groupGuildFunctionHeatmapTopAsvs = groupGuildFunctionConfig.heatmap_top_asvs != null ? (groupGuildFunctionConfig.heatmap_top_asvs as int) : 80
    def groupGuildFunctionPcaLabelMinMeanAbundancePct = groupGuildFunctionConfig.pca_label_min_mean_abundance_pct != null ? (groupGuildFunctionConfig.pca_label_min_mean_abundance_pct as double) : 1.0
    if( !Double.isFinite(groupGuildFunctionPcaLabelMinMeanAbundancePct) || groupGuildFunctionPcaLabelMinMeanAbundancePct < 0 || groupGuildFunctionPcaLabelMinMeanAbundancePct >= 100 ) {
        exit 1, "group_guild_function.pca_label_min_mean_abundance_pct must be >= 0 and < 100"
    }
    def groupGuildFunctionPcaScores = resolveOptionalPath(
        groupGuildFunctionConfig.pca_scores, configRoot
    )
    def groupGuildFunctionPcaExplained = resolveOptionalPath(
        groupGuildFunctionConfig.pca_explained, configRoot
    )
    def groupGuildFunctionHybridAssignments = resolveOptionalPath(
        groupGuildFunctionConfig.hybrid_assignments, configRoot
    )
    def groupGuildFunctionHybridCentroids = resolveOptionalPath(
        groupGuildFunctionConfig.hybrid_centroids, configRoot
    )
    def groupGuildFunctionPcaOverlayPaths = [
        groupGuildFunctionPcaScores,
        groupGuildFunctionPcaExplained,
        groupGuildFunctionHybridAssignments,
        groupGuildFunctionHybridCentroids,
    ]
    boolean groupGuildFunctionPcaOverlayEnabled =
        groupGuildFunctionPcaOverlayPaths.every { it }
    if( groupGuildFunctionPcaOverlayPaths.any { it } && !groupGuildFunctionPcaOverlayEnabled ) {
        exit 1, "group_guild_function PCA overlay requires pca_scores, pca_explained, hybrid_assignments, and hybrid_centroids"
    }
    if( groupGuildFunctionPcaOverlayEnabled ) {
        groupGuildFunctionPcaOverlayPaths.each { requiredPath ->
            if( !file(requiredPath).exists() ) {
                exit 1, "Group-guild-function PCA overlay input not found: ${requiredPath}"
            }
        }
    }

    def powerAnalysisConfig = config.group_power_analysis ?: (config.power_analysis ?: [:])
    boolean powerAnalysisRequested = powerAnalysisConfig.containsKey('enabled') ? (powerAnalysisConfig.enabled as boolean) : false
    if( powerAnalysisRequested && !metadataPlotsEnabled ) {
        exit 1, "power_analysis.enabled requires metadata_plots.enabled to be true"
    }
    boolean powerAnalysisEnabled = powerAnalysisRequested
    def powerAnalysisOutputDir = powerAnalysisConfig.output_dir ?: 'power_analysis'
    def powerAnalysisOutputDirAbs = resolveOutputRelative(powerAnalysisOutputDir.toString(), outputDir)
    def powerAnalysisSampleCol = powerAnalysisConfig.sample_col ?: metadataPlotsSampleCol
    def powerAnalysisPatientCol = powerAnalysisConfig.patient_col ? powerAnalysisConfig.patient_col.toString().trim() : 'Participant_ID'
    def powerAnalysisCaseCol = powerAnalysisConfig.case_col ? powerAnalysisConfig.case_col.toString().trim() : 'Case'
    def powerAnalysisTypeCol = powerAnalysisConfig.type_col ?: metadataPlotsTypeCol
    def powerAnalysisSampleSizesCancerRaw = powerAnalysisConfig.sample_sizes_cancer ?: '6,8,10,15,20,25,30'
    def powerAnalysisSampleSizesCancer = powerAnalysisSampleSizesCancerRaw instanceof List ?
        powerAnalysisSampleSizesCancerRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        powerAnalysisSampleSizesCancerRaw.toString().trim()
    def powerAnalysisSampleSizesStypeRaw = powerAnalysisConfig.sample_sizes_stype ?: '10,15,20,25,30,40,50'
    def powerAnalysisSampleSizesStype = powerAnalysisSampleSizesStypeRaw instanceof List ?
        powerAnalysisSampleSizesStypeRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        powerAnalysisSampleSizesStypeRaw.toString().trim()
    def powerAnalysisNSimulations = powerAnalysisConfig.n_simulations ? (powerAnalysisConfig.n_simulations as int) : 1000
    def powerAnalysisNPerm = powerAnalysisConfig.n_perm ? (powerAnalysisConfig.n_perm as int) : 199
    def powerAnalysisAlpha = powerAnalysisConfig.alpha != null ? (powerAnalysisConfig.alpha as double) : 0.05d
    def powerAnalysisSeed = powerAnalysisConfig.seed ? (powerAnalysisConfig.seed as int) : 42
    boolean powerAnalysisSkipEstimate = powerAnalysisConfig.containsKey('skip_estimate') ? (powerAnalysisConfig.skip_estimate as boolean) : false
    boolean powerAnalysisSkipPlot = powerAnalysisConfig.containsKey('skip_plot') ? (powerAnalysisConfig.skip_plot as boolean) : false
    def powerAnalysisTransform = powerAnalysisConfig.transform ? powerAnalysisConfig.transform.toString().trim().toLowerCase() : 'none'
    if( !['none', 'rclr'].contains(powerAnalysisTransform) ) {
        exit 1, "power_analysis.transform must be one of: none, rclr"
    }
    boolean powerAnalysisKeepContralateralInCancer = powerAnalysisConfig.containsKey('keep_contralateral_in_cancer') ? (powerAnalysisConfig.keep_contralateral_in_cancer as boolean) : false
    def powerAnalysisContralateralTypesRaw = powerAnalysisConfig.contralateral_sample_types ?: 'Lung Brush,BAL'
    def powerAnalysisContralateralTypes = powerAnalysisContralateralTypesRaw instanceof List ?
        powerAnalysisContralateralTypesRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        powerAnalysisContralateralTypesRaw.toString().trim()
    def powerAnalysisIndicspeciesDir = powerAnalysisConfig.indicspecies_dir ? resolveOptionalPath(powerAnalysisConfig.indicspecies_dir, configRoot) : indicspeciesOutputDirAbs

    def taxonomyPatientAwareConfig = config.taxonomy_group_association ?: (config.taxonomy_patient_aware ?: [:])
    boolean taxonomyPatientAwareRequested = taxonomyPatientAwareConfig.containsKey('enabled') ? (taxonomyPatientAwareConfig.enabled as boolean) : false
    if( taxonomyPatientAwareRequested && !metadataPlotsEnabled ) {
        exit 1, "taxonomy_patient_aware.enabled requires metadata_plots.enabled to be true"
    }
    boolean taxonomyPatientAwareEnabled = taxonomyPatientAwareRequested
    def taxonomyPatientAwareOutputDir = taxonomyPatientAwareConfig.output_dir ?: (config.taxonomy_group_association ? 'taxonomy_group_association' : 'taxonomy_patient_aware')
    def taxonomyPatientAwareOutputDirAbs = resolveOutputRelative(taxonomyPatientAwareOutputDir.toString(), outputDir)
    def taxonomyPatientAwareSampleCol = taxonomyPatientAwareConfig.sample_col ?: metadataPlotsSampleCol
    def taxonomyPatientAwarePatientCol = taxonomyPatientAwareConfig.subject_col ? taxonomyPatientAwareConfig.subject_col.toString().trim() : (taxonomyPatientAwareConfig.patient_col ? taxonomyPatientAwareConfig.patient_col.toString().trim() : 'Participant_ID')
    def taxonomyPatientAwareCaseCol = taxonomyPatientAwareConfig.comparison_col ? taxonomyPatientAwareConfig.comparison_col.toString().trim() : (taxonomyPatientAwareConfig.case_col ? taxonomyPatientAwareConfig.case_col.toString().trim() : 'Case')
    def taxonomyPatientAwareTypeCol = taxonomyPatientAwareConfig.group_col ?: (taxonomyPatientAwareConfig.type_col ?: metadataPlotsTypeCol)
    def taxonomyPatientAwareCountCol = taxonomyPatientAwareConfig.count_col ? taxonomyPatientAwareConfig.count_col.toString().trim() : 'count'
    def taxonomyPatientAwareTaxLevelsRaw = taxonomyPatientAwareConfig.tax_levels ?: 'Phylum,Family'
    def taxonomyPatientAwareTaxLevels = taxonomyPatientAwareTaxLevelsRaw instanceof List ?
        taxonomyPatientAwareTaxLevelsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        taxonomyPatientAwareTaxLevelsRaw.toString().trim()
    def taxonomyPatientAwareSampleTypesRaw = taxonomyPatientAwareConfig.containsKey('sample_types') ? taxonomyPatientAwareConfig.sample_types : 'Oral Rinse,BAL,Lung Brush'
    def taxonomyPatientAwareSampleTypes = taxonomyPatientAwareSampleTypesRaw instanceof List ?
        taxonomyPatientAwareSampleTypesRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        taxonomyPatientAwareSampleTypesRaw.toString().trim()
    def taxonomyPatientAwareComparisonGroupsRaw = taxonomyPatientAwareConfig.comparison_groups ?: ''
    def taxonomyPatientAwareComparisonGroups = taxonomyPatientAwareComparisonGroupsRaw instanceof List ?
        taxonomyPatientAwareComparisonGroupsRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        taxonomyPatientAwareComparisonGroupsRaw.toString().trim()
    boolean taxonomyPatientAwareRunComparison = taxonomyPatientAwareConfig.containsKey('run_comparison') ?
        (taxonomyPatientAwareConfig.run_comparison as boolean) :
        (!config.taxonomy_group_association || taxonomyPatientAwareComparisonGroups.toString().trim().length() > 0)
    def taxonomyPatientAwareMinPrevalence = taxonomyPatientAwareConfig.min_prevalence != null ? (taxonomyPatientAwareConfig.min_prevalence as double) : 0.10d
    boolean taxonomyPatientAwareExcludeContralateral = taxonomyPatientAwareConfig.containsKey('exclude_contralateral_in_cancer') ? (taxonomyPatientAwareConfig.exclude_contralateral_in_cancer as boolean) : true
    def taxonomyPatientAwareContralateralCol = taxonomyPatientAwareConfig.contralateral_col ? taxonomyPatientAwareConfig.contralateral_col.toString().trim() : 'lung_status'
    def taxonomyPatientAwareCancerSiteCol = taxonomyPatientAwareConfig.cancer_site_col ? taxonomyPatientAwareConfig.cancer_site_col.toString().trim() : 'Cancer_Site'
    def taxonomyPatientAwareLungSideCol = taxonomyPatientAwareConfig.lung_side_col ? taxonomyPatientAwareConfig.lung_side_col.toString().trim() : 'lung_code'
    def taxonomyPatientAwareContralateralValue = taxonomyPatientAwareConfig.contralateral_value ? taxonomyPatientAwareConfig.contralateral_value.toString().trim() : 'Contralateral'
    def taxonomyPatientAwareContralateralTypesRaw = taxonomyPatientAwareConfig.contralateral_sample_types ?: 'Lung Brush,BAL'
    def taxonomyPatientAwareContralateralTypes = taxonomyPatientAwareContralateralTypesRaw instanceof List ?
        taxonomyPatientAwareContralateralTypesRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        taxonomyPatientAwareContralateralTypesRaw.toString().trim()
    boolean taxonomyPatientAwareSkipOmnibus = taxonomyPatientAwareConfig.containsKey('skip_omnibus') ? (taxonomyPatientAwareConfig.skip_omnibus as boolean) : false
    def taxonomyPatientAwareTransform = taxonomyPatientAwareConfig.transform ? taxonomyPatientAwareConfig.transform.toString().trim().toLowerCase() : 'none'
    if( !['none', 'rclr'].contains(taxonomyPatientAwareTransform) ) {
        exit 1, "taxonomy_patient_aware.transform must be one of: none, rclr"
    }
    def taxonomyPatientAwareAlpha = taxonomyPatientAwareConfig.alpha != null ? (taxonomyPatientAwareConfig.alpha as double) : 0.05d
    def taxonomyPatientAwareTopN = taxonomyPatientAwareConfig.top_n ? (taxonomyPatientAwareConfig.top_n as int) : 12
    def taxonomyPatientAwareTypePalette = taxonomyPatientAwareConfig.type_palette ?: (indicspeciesGroupPaletteMap[taxonomyPatientAwareTypeCol] ?: '')
    def taxonomyPatientAwareCasePalette = taxonomyPatientAwareConfig.case_palette ?: (indicspeciesGroupPaletteMap[taxonomyPatientAwareCaseCol] ?: '')

    def lungStatusAnalysisConfig = config.paired_group_contrast ?: (config.lung_status_analysis ?: [:])
    boolean lungStatusAnalysisRequested = lungStatusAnalysisConfig.containsKey('enabled') ? (lungStatusAnalysisConfig.enabled as boolean) : false
    if( lungStatusAnalysisRequested && !metadataPlotsEnabled ) {
        exit 1, "lung_status_analysis.enabled requires metadata_plots.enabled to be true"
    }
    boolean lungStatusAnalysisEnabled = lungStatusAnalysisRequested
    def lungStatusAnalysisOutputDir = lungStatusAnalysisConfig.output_dir ?: (config.paired_group_contrast ? 'paired_group_contrast' : 'lung_status_analysis')
    def lungStatusAnalysisOutputDirAbs = resolveOutputRelative(lungStatusAnalysisOutputDir.toString(), outputDir)
    def lungStatusAnalysisSampleCol = lungStatusAnalysisConfig.sample_col ?: metadataPlotsSampleCol
    def lungStatusAnalysisTypeCol = lungStatusAnalysisConfig.type_col ?: metadataPlotsTypeCol
    def lungStatusAnalysisSampleTypesRaw = lungStatusAnalysisConfig.containsKey('sample_types') ? lungStatusAnalysisConfig.sample_types : 'Lung Brush,BAL'
    def lungStatusAnalysisSampleTypes = lungStatusAnalysisSampleTypesRaw instanceof List ?
        lungStatusAnalysisSampleTypesRaw.collect { it.toString().trim() }.findAll { it }.join(',') :
        lungStatusAnalysisSampleTypesRaw.toString().trim()
    def lungStatusAnalysisCaseCol = lungStatusAnalysisConfig.case_col ? lungStatusAnalysisConfig.case_col.toString().trim() : 'Case'
    def lungStatusAnalysisPatientCol = lungStatusAnalysisConfig.subject_col ? lungStatusAnalysisConfig.subject_col.toString().trim() : (lungStatusAnalysisConfig.patient_col ? lungStatusAnalysisConfig.patient_col.toString().trim() : 'Participant_ID')
    def lungStatusAnalysisCancerSiteCol = lungStatusAnalysisConfig.cancer_site_col ? lungStatusAnalysisConfig.cancer_site_col.toString().trim() : 'Cancer_Site'
    def lungStatusAnalysisLungCodeCol = lungStatusAnalysisConfig.lung_code_col ? lungStatusAnalysisConfig.lung_code_col.toString().trim() : 'lung_code'
    def lungStatusAnalysisTumorSideCol = lungStatusAnalysisConfig.tumor_side_col ? lungStatusAnalysisConfig.tumor_side_col.toString().trim() : 'TumorSide'
    def lungStatusAnalysisContralateralCol = lungStatusAnalysisConfig.contralateral_col ? lungStatusAnalysisConfig.contralateral_col.toString().trim() : 'Contralateral'
    def lungStatusAnalysisHealthyCol = lungStatusAnalysisConfig.healthy_col ? lungStatusAnalysisConfig.healthy_col.toString().trim() : 'Healthy'
    def lungStatusAnalysisStatusCol = lungStatusAnalysisConfig.lung_status_col ? lungStatusAnalysisConfig.lung_status_col.toString().trim() : 'lung_status'
    def lungStatusAnalysisStatusAValue = lungStatusAnalysisConfig.status_a_value ? lungStatusAnalysisConfig.status_a_value.toString().trim() : 'TumorSide'
    def lungStatusAnalysisStatusBValue = lungStatusAnalysisConfig.status_b_value ? lungStatusAnalysisConfig.status_b_value.toString().trim() : 'Contralateral'
    def lungStatusAnalysisReferenceStatusValue = lungStatusAnalysisConfig.reference_status_value ? lungStatusAnalysisConfig.reference_status_value.toString().trim() : 'Healthy'
    def lungStatusAnalysisPermutations = lungStatusAnalysisConfig.permutations ? (lungStatusAnalysisConfig.permutations as int) : 9999
    def lungStatusAnalysisSeed = lungStatusAnalysisConfig.seed ? (lungStatusAnalysisConfig.seed as int) : 1
    return [
        metadataPlotsEnabled: metadataPlotsEnabled,
        metadataPlotsMetadataPath: metadataPlotsMetadataPath,
        metadataPlotsSampleCol: metadataPlotsSampleCol,
        metadataPlotsTypeCol: metadataPlotsTypeCol,
        metadataPlotsColorCol: metadataPlotsColorCol,
        metadataPlotsGroupOrder: metadataPlotsGroupOrder,
        indicspeciesSampleCol: indicspeciesSampleCol,
        indicspeciesPerms: indicspeciesPerms,
        indicspeciesSeed: indicspeciesSeed,
        indicspeciesQThreshold: indicspeciesQThreshold,
        indicspeciesMinN: indicspeciesMinN,
        indicspeciesBlockCol: indicspeciesBlockCol,
        indicspeciesBlockedCols: indicspeciesBlockedCols,
        indicspeciesStratifiedSpecsArg: indicspeciesStratifiedSpecsArg,
        indicspeciesGroup1: indicspeciesGroup1,
        indicspeciesGroup2: indicspeciesGroup2,
        indicspeciesOutputDirAbs: indicspeciesOutputDirAbs,
        indicspeciesPlotPairsMode: indicspeciesPlotPairsMode,
        indicspeciesPlotOutputDirAbs: indicspeciesPlotOutputDirAbs,
        indicspeciesPlotVennPath: indicspeciesPlotVennPath,
        indicspeciesPlotTaxonomyPath: indicspeciesPlotTaxonomyPath,
        indicspeciesColorCol: indicspeciesColorCol,
        indicspeciesGroupPaletteJson: indicspeciesGroupPaletteJson,
        indicspeciesGroupOrderJson: indicspeciesGroupOrderJson,
        indicspeciesFocusLabelJson: indicspeciesFocusLabelJson,
        indicspeciesAlignedOutputDirAbs: indicspeciesAlignedOutputDirAbs,
        indicspeciesAlignedAlpha: indicspeciesAlignedAlpha,
        indicspeciesAlignedMinStat: indicspeciesAlignedMinStat,
        indicspeciesAlignedTopN: indicspeciesAlignedTopN,
        vocCorrelationEnabled: vocCorrelationEnabled,
        vocCorrelationVocTablePath: vocCorrelationVocTablePath,
        vocCorrelationOutputDirAbs: vocCorrelationOutputDirAbs,
        vocCorrelationVocSampleCol: vocCorrelationVocSampleCol,
        vocCorrelationSampleIdMode: vocCorrelationSampleIdMode,
        vocCorrelationMetadataSampleCol: vocCorrelationMetadataSampleCol,
        vocCorrelationTypeCol: vocCorrelationTypeCol,
        vocCorrelationPatientCol: vocCorrelationPatientCol,
        vocCorrelationCaseCol: vocCorrelationCaseCol,
        vocCorrelationSampleTypes: vocCorrelationSampleTypes,
        vocCorrelationUseLegacySubset: vocCorrelationUseLegacySubset,
        vocCorrelationVocCols: vocCorrelationVocCols,
        vocCorrelationDirection: vocCorrelationDirection,
        vocCorrelationCasePalette: vocCorrelationCasePalette,
        vocCorrelationIsaPalette: vocCorrelationIsaPalette,
        measurementAssociationEnabled: measurementAssociationEnabled,
        measurementAssociationOutputDirAbs: measurementAssociationOutputDirAbs,
        measurementAssociationTablePath: measurementAssociationTablePath,
        measurementAssociationSampleCol: measurementAssociationSampleCol,
        measurementAssociationAsvIdCol: measurementAssociationAsvIdCol,
        measurementAssociationMeasurementSampleCol: measurementAssociationMeasurementSampleCol,
        measurementAssociationMetadataJoinCols: measurementAssociationMetadataJoinCols,
        measurementAssociationMeasurementJoinCols: measurementAssociationMeasurementJoinCols,
        measurementAssociationAliasesJson: measurementAssociationAliasesJson,
        measurementAssociationCols: measurementAssociationCols,
        measurementAssociationExcludeCols: measurementAssociationExcludeCols,
        measurementAssociationGroupCol: measurementAssociationGroupCol,
        measurementAssociationGroupPalette: measurementAssociationGroupPalette,
        measurementAssociationSubsetSource: measurementAssociationSubsetSource,
        measurementAssociationMaxAsvs: measurementAssociationMaxAsvs,
        measurementAssociationMinTotal: measurementAssociationMinTotal,
        measurementAssociationMinPrevalence: measurementAssociationMinPrevalence,
        measurementAssociationTopCorrelations: measurementAssociationTopCorrelations,
        measurementAssociationDirection: measurementAssociationDirection,
        measurementAssociationMethods: measurementAssociationMethods,
        measurementAssociationOrdinationCols: measurementAssociationOrdinationCols,
        measurementAssociationPermutations: measurementAssociationPermutations,
        measurementAssociationTopVectors: measurementAssociationTopVectors,
        measurementAssociationFormats: measurementAssociationFormats,
        moduleMeasurementAssociationEnabled: moduleMeasurementAssociationEnabled,
        moduleMeasurementSparseCols: moduleMeasurementSparseCols,
        moduleMeasurementCruiseCol: moduleMeasurementCruiseCol,
        moduleMeasurementDepthCol: moduleMeasurementDepthCol,
        moduleMeasurementQThreshold: moduleMeasurementQThreshold,
        moduleMeasurementNetworkMetricTopN: moduleMeasurementNetworkMetricTopN,
        moduleMeasurementPcaScores: moduleMeasurementPcaScores,
        moduleMeasurementPcaLoadings: moduleMeasurementPcaLoadings,
        moduleMeasurementPcaExplained: moduleMeasurementPcaExplained,
        moduleMeasurementHybridAssignments: moduleMeasurementHybridAssignments,
        moduleMeasurementHybridCentroids: moduleMeasurementHybridCentroids,
        titanEnabled: titanEnabled,
        titanOutputDirAbs: titanOutputDirAbs,
        titanVariables: titanVariables,
        titanTranspose: titanTranspose,
        titanMinSplit: titanMinSplit,
        titanMinimumSamples: titanMinimumSamples,
        titanMinimumOccurrence: titanMinimumOccurrence,
        titanMinimumPrevalence: titanMinimumPrevalence,
        titanMinimumMeanRelativeAbundance: titanMinimumMeanRelativeAbundance,
        titanPermutations: titanPermutations,
        titanBootstrapCount: titanBootstrapCount,
        titanSeed: titanSeed,
        titanSourceRepository: titanSourceRepository,
        titanSourceRevision: titanSourceRevision,
        titanExpectedVersion: titanExpectedVersion,
        titanImax: titanImax,
        titanIvTotal: titanIvTotal,
        titanPurityCutoff: titanPurityCutoff,
        titanReliabilityCutoff: titanReliabilityCutoff,
        titanNcpus: titanNcpus,
        titanMemory: titanMemory,
        titanPlotEnabled: titanPlotEnabled,
        titanRankedTopN: titanRankedTopN,
        titanFormats: titanFormats,
        microbialCompartmentEnabled: microbialCompartmentEnabled,
        microbialCompartmentOutputDirAbs: microbialCompartmentOutputDirAbs,
        microbialCompartmentSampleCol: microbialCompartmentSampleCol,
        microbialCompartmentTranspose: microbialCompartmentTranspose,
        microbialCompartmentMinimumPrevalence: microbialCompartmentMinimumPrevalence,
        microbialCompartmentMinimumMaxRelativeAbundance: microbialCompartmentMinimumMaxRelativeAbundance,
        microbialCompartmentZeroReplacement: microbialCompartmentZeroReplacement,
        microbialCompartmentMultiplicativeDelta: microbialCompartmentMultiplicativeDelta,
        microbialCompartmentPseudocount: microbialCompartmentPseudocount,
        microbialCompartmentKMin: microbialCompartmentKMin,
        microbialCompartmentKMax: microbialCompartmentKMax,
        microbialCompartmentMinClusterSize: microbialCompartmentMinClusterSize,
        microbialCompartmentMinClusterFraction: microbialCompartmentMinClusterFraction,
        microbialCompartmentMinMeanSilhouette: microbialCompartmentMinMeanSilhouette,
        microbialCompartmentMinStabilityAri: microbialCompartmentMinStabilityAri,
        microbialCompartmentMinClusterJaccard: microbialCompartmentMinClusterJaccard,
        microbialCompartmentStabilityReplicates: microbialCompartmentStabilityReplicates,
        microbialCompartmentStabilitySampleFraction: microbialCompartmentStabilitySampleFraction,
        microbialCompartmentStabilityBlockCol: microbialCompartmentStabilityBlockCol,
        microbialCompartmentStabilityStratumCol: microbialCompartmentStabilityStratumCol,
        microbialCompartmentStabilityPrimaryQuantile: microbialCompartmentStabilityPrimaryQuantile,
        microbialCompartmentStabilityNearTieTolerance: microbialCompartmentStabilityNearTieTolerance,
        microbialCompartmentBlockedRobustnessEnabled: microbialCompartmentBlockedRobustnessEnabled,
        microbialCompartmentPredictionStrengthEnabled: microbialCompartmentPredictionStrengthEnabled,
        microbialCompartmentHierarchicalEnabled: microbialCompartmentHierarchicalEnabled,
        microbialCompartmentSeed: microbialCompartmentSeed,
        microbialCompartmentNcpus: microbialCompartmentNcpus,
        microbialCompartmentEnvironmentalCols: microbialCompartmentEnvironmentalCols,
        microbialCompartmentDepthCol: microbialCompartmentDepthCol,
        microbialCompartmentDateCol: microbialCompartmentDateCol,
        microbialCompartmentDominantTopN: microbialCompartmentDominantTopN,
        microbialCompartmentPlotEnabled: microbialCompartmentPlotEnabled,
        microbialCompartmentFormats: microbialCompartmentFormats,
        indicspeciesRequiresMicrobialCompartments: indicspeciesRequiresMicrobialCompartments,
        indicspeciesRunEarly: indicspeciesRunEarly,
        microbialStateInterpretationEnabled: microbialStateInterpretationEnabled,
        microbialStateInterpretationOutputDirAbs: microbialStateInterpretationOutputDirAbs,
        microbialStateInterpretationSampleCol: microbialStateInterpretationSampleCol,
        microbialStateInterpretationCruiseCol: microbialStateInterpretationCruiseCol,
        microbialStateInterpretationDepthCol: microbialStateInterpretationDepthCol,
        microbialStateInterpretationSeasonCol: microbialStateInterpretationSeasonCol,
        microbialStateInterpretationYearCol: microbialStateInterpretationYearCol,
        microbialStateInterpretationDateCol: microbialStateInterpretationDateCol,
        microbialStateInterpretationRenewalCol: microbialStateInterpretationRenewalCol,
        microbialStateInterpretationRenewalLevels: microbialStateInterpretationRenewalLevels,
        microbialStateInterpretationEnvironmentalCols: microbialStateInterpretationEnvironmentalCols,
        microbialStateInterpretationHybridCol: microbialStateInterpretationHybridCol,
        microbialStateInterpretationPermutations: microbialStateInterpretationPermutations,
        microbialStateInterpretationBootstrapReplicates: microbialStateInterpretationBootstrapReplicates,
        microbialStateInterpretationSeed: microbialStateInterpretationSeed,
        microbialStateInterpretationHybridPalette: microbialStateInterpretationHybridPalette,
        microbialStateInterpretationHybridOrder: microbialStateInterpretationHybridOrder,
        microbialStateInterpretationMcPalette: microbialStateInterpretationMcPalette,
        microbialStateInterpretationDepthPalette: microbialStateInterpretationDepthPalette,
        microbialStateInterpretationMcOrder: microbialStateInterpretationMcOrder,
        microbialStateInterpretationAsvTopN: microbialStateInterpretationAsvTopN,
        microbialStateInterpretationSelectedAsvs: microbialStateInterpretationSelectedAsvs,
        microbialStateInterpretationAsvSelectionMetric: microbialStateInterpretationAsvSelectionMetric,
        microbialStateInterpretationMaximumDepth: microbialStateInterpretationMaximumDepth,
        microbialStateInterpretationDepthStep: microbialStateInterpretationDepthStep,
        microbialStateInterpretationTimeSubdivisions: microbialStateInterpretationTimeSubdivisions,
        microbialStateInterpretationTimeSigma: microbialStateInterpretationTimeSigma,
        microbialStateInterpretationDepthSigma: microbialStateInterpretationDepthSigma,
        microbialStateInterpretationExcludeCurtainPattern: microbialStateInterpretationExcludeCurtainPattern,
        microbialStateInterpretationRenewalEvents: microbialStateInterpretationRenewalEvents,
        microbialStateInterpretationRenewalDateCol: microbialStateInterpretationRenewalDateCol,
        microbialStateInterpretationFormats: microbialStateInterpretationFormats,
        ecologicalContextAtlasEnabled: ecologicalContextAtlasEnabled,
        ecologicalContextAtlasOutputDirAbs: ecologicalContextAtlasOutputDirAbs,
        ecologicalContextAtlasSampleCol: ecologicalContextAtlasSampleCol,
        ecologicalContextAtlasCruiseCol: ecologicalContextAtlasCruiseCol,
        ecologicalContextAtlasDepthCol: ecologicalContextAtlasDepthCol,
        ecologicalContextAtlasDateCol: ecologicalContextAtlasDateCol,
        ecologicalContextAtlasContextCols: ecologicalContextAtlasContextCols,
        ecologicalContextAtlasNetworkVariables: ecologicalContextAtlasNetworkVariables,
        ecologicalContextAtlasLinkageQThreshold: ecologicalContextAtlasLinkageQThreshold,
        ecologicalContextAtlasSelectedAsvs: ecologicalContextAtlasSelectedAsvs,
        ecologicalContextAtlasContextSheetSource: ecologicalContextAtlasContextSheetSource,
        ecologicalContextAtlasRnaDnaLog2TpmRatio: ecologicalContextAtlasRnaDnaLog2TpmRatio,
        ecologicalContextAtlasRnaDnaLog2TpmSe: ecologicalContextAtlasRnaDnaLog2TpmSe,
        ecologicalContextAtlasRnaDnaMaxLog2Se: ecologicalContextAtlasRnaDnaMaxLog2Se,
        ecologicalContextAtlasRnaDnaMagIdMode: ecologicalContextAtlasRnaDnaMagIdMode,
        ecologicalContextAtlasFormats: ecologicalContextAtlasFormats,
        communityTurnoverEnabled: communityTurnoverEnabled,
        communityTurnoverOutputDirAbs: communityTurnoverOutputDirAbs,
        communityTurnoverSampleCol: communityTurnoverSampleCol,
        communityTurnoverProfileCol: communityTurnoverProfileCol,
        communityTurnoverDateCol: communityTurnoverDateCol,
        communityTurnoverDepthCol: communityTurnoverDepthCol,
        communityTurnoverTranspose: communityTurnoverTranspose,
        communityTurnoverMinimumPrevalence: communityTurnoverMinimumPrevalence,
        communityTurnoverMinimumMaxRelativeAbundance: communityTurnoverMinimumMaxRelativeAbundance,
        communityTurnoverZeroReplacement: communityTurnoverZeroReplacement,
        communityTurnoverMultiplicativeDelta: communityTurnoverMultiplicativeDelta,
        communityTurnoverPseudocount: communityTurnoverPseudocount,
        communityTurnoverDistanceMetrics: communityTurnoverDistanceMetrics,
        communityTurnoverPrimaryMetric: communityTurnoverPrimaryMetric,
        communityTurnoverFixedDepthMinProfileFraction: communityTurnoverFixedDepthMinProfileFraction,
        communityTurnoverMinimumTimeDifferenceDays: communityTurnoverMinimumTimeDifferenceDays,
        communityTurnoverEnvironmentalCols: communityTurnoverEnvironmentalCols,
        communityTurnoverPrimaryEnvironmentalCol: communityTurnoverPrimaryEnvironmentalCol,
        communityTurnoverLcbdPermutations: communityTurnoverLcbdPermutations,
        communityTurnoverBoundaryPermutations: communityTurnoverBoundaryPermutations,
        communityTurnoverBoundaryBootstrapReplicates: communityTurnoverBoundaryBootstrapReplicates,
        communityTurnoverWithinProfileLcbdEnabled: communityTurnoverWithinProfileLcbdEnabled,
        communityTurnoverWithinProfileLcbdMinSamples: communityTurnoverWithinProfileLcbdMinSamples,
        communityTurnoverSeed: communityTurnoverSeed,
        communityTurnoverPlotEnabled: communityTurnoverPlotEnabled,
        communityTurnoverFormats: communityTurnoverFormats,
        communityPredictorEnabled: communityPredictorEnabled,
        communityPredictorOutputDirAbs: communityPredictorOutputDirAbs,
        communityPredictorSampleCol: communityPredictorSampleCol,
        communityPredictorCruiseCol: communityPredictorCruiseCol,
        communityPredictorYearCol: communityPredictorYearCol,
        communityPredictorDateCol: communityPredictorDateCol,
        communityPredictorSeasonCol: communityPredictorSeasonCol,
        communityPredictorDepthCol: communityPredictorDepthCol,
        communityPredictorCruiseDepthMinPrevalence: communityPredictorCruiseDepthMinPrevalence,
        communityPredictorPeaCol: communityPredictorPeaCol,
        communityPredictorCentroidCol: communityPredictorCentroidCol,
        communityPredictorCruiseGroupCol: communityPredictorCruiseGroupCol,
        communityPredictorCruiseGroupProbabilityCol: communityPredictorCruiseGroupProbabilityCol,
        communityPredictorCruiseGroupUncertainCol: communityPredictorCruiseGroupUncertainCol,
        communityPredictorRenewalGroupCol: communityPredictorRenewalGroupCol,
        communityPredictorO2Col: communityPredictorO2Col,
        communityPredictorGmmCol: communityPredictorGmmCol,
        communityPredictorHybridCol: communityPredictorHybridCol,
        communityPredictorSourceCol: communityPredictorSourceCol,
        communityPredictorObservedLabel: communityPredictorObservedLabel,
        communityPredictorMinGroupN: communityPredictorMinGroupN,
        communityPredictorPermutations: communityPredictorPermutations,
        communityPredictorSeed: communityPredictorSeed,
        communityPredictorFormats: communityPredictorFormats,
        groupingDiagnosticsEnabled: groupingDiagnosticsEnabled,
        groupingDiagnosticsRequiresMicrobialCompartments: groupingDiagnosticsRequiresMicrobialCompartments,
        groupingDiagnosticsRunEarly: groupingDiagnosticsRunEarly,
        groupingDiagnosticsOutputDirAbs: groupingDiagnosticsOutputDirAbs,
        groupingDiagnosticsSampleCol: groupingDiagnosticsSampleCol,
        groupingDiagnosticsGroupCols: groupingDiagnosticsGroupCols,
        groupingDiagnosticsCruiseGroupCols: groupingDiagnosticsCruiseGroupCols,
        groupingDiagnosticsCruiseCol: groupingDiagnosticsCruiseCol,
        groupingDiagnosticsDepthCol: groupingDiagnosticsDepthCol,
        groupingDiagnosticsCruiseDepthMinPrevalence: groupingDiagnosticsCruiseDepthMinPrevalence,
        groupingDiagnosticsBaselineGroup: groupingDiagnosticsBaselineGroup,
        groupingDiagnosticsPrimaryGroup: groupingDiagnosticsPrimaryGroup,
        groupingDiagnosticsPaletteJson: groupingDiagnosticsPaletteJson,
        groupingDiagnosticsOrderJson: groupingDiagnosticsOrderJson,
        groupingDiagnosticsMetrics: groupingDiagnosticsMetrics,
        groupingDiagnosticsTransform: groupingDiagnosticsTransform,
        groupingDiagnosticsPermutations: groupingDiagnosticsPermutations,
        groupingDiagnosticsRandomState: groupingDiagnosticsRandomState,
        groupingDiagnosticsFormats: groupingDiagnosticsFormats,
        groupingDiagnosticsSoftLabelEnabled: groupingDiagnosticsSoftLabelEnabled,
        groupingDiagnosticsSoftLabelK: groupingDiagnosticsSoftLabelK,
        groupingDiagnosticsSoftLabelTargetCols: groupingDiagnosticsSoftLabelTargetCols,
        groupingDiagnosticsSoftLabelExcludeLabels: groupingDiagnosticsSoftLabelExcludeLabels,
        groupingDiagnosticsSoftLabelMinClassSamples: groupingDiagnosticsSoftLabelMinClassSamples,
        groupingDiagnosticsSoftLabelDistanceQuantile: groupingDiagnosticsSoftLabelDistanceQuantile,
        groupingDiagnosticsApplySoftLabels: groupingDiagnosticsApplySoftLabels,
        groupingDiagnosticsSoftLabelTargetCol: groupingDiagnosticsSoftLabelTargetCol,
        groupingDiagnosticsSoftLabelMinConfidence: groupingDiagnosticsSoftLabelMinConfidence,
        groupingDiagnosticsSoftLabelMinNeighborAgreement: groupingDiagnosticsSoftLabelMinNeighborAgreement,
        groupingDiagnosticsSoftLabelMinCvBalancedAccuracy: groupingDiagnosticsSoftLabelMinCvBalancedAccuracy,
        groupingDiagnosticsPowerEnabled: groupingDiagnosticsPowerEnabled,
        groupingDiagnosticsPowerSizes: groupingDiagnosticsPowerSizes,
        groupingDiagnosticsPowerSimulations: groupingDiagnosticsPowerSimulations,
        groupingDiagnosticsPowerPermutations: groupingDiagnosticsPowerPermutations,
        groupingDiagnosticsPowerAlpha: groupingDiagnosticsPowerAlpha,
        groupingDiagnosticsPowerMinGroups: groupingDiagnosticsPowerMinGroups,
        clustermapsOutputDirAbs: clustermapsOutputDirAbs,
        clustermapsMitoOutputDirAbs: clustermapsMitoOutputDirAbs,
        clustermapsMitoInputPath: clustermapsMitoInputPath,
        clustermapsSampleCol: clustermapsSampleCol,
        clustermapsSampleCodeCol: clustermapsSampleCodeCol,
        clustermapsAsvIdCol: clustermapsAsvIdCol,
        clustermapsGroup1Col: clustermapsGroup1Col,
        clustermapsGroup2Col: clustermapsGroup2Col,
        clustermapsGroup3Col: clustermapsGroup3Col,
        clustermapsExcludeGroup1: clustermapsExcludeGroup1,
        clustermapsGroup1Palette: clustermapsGroup1Palette,
        clustermapsGroup2Palette: clustermapsGroup2Palette,
        clustermapsGroup3Palette: clustermapsGroup3Palette,
        clustermapsRanks: clustermapsRanks,
        clustermapsTopN: clustermapsTopN,
        clustermapsPlotAsvLevel: clustermapsPlotAsvLevel,
        clustermapsCountCol: clustermapsCountCol,
        clustermapsIsaMinStat: clustermapsIsaMinStat,
        clustermapsIsaSignificanceCols: clustermapsIsaSignificanceCols,
        clustermapsIsaStatCols: clustermapsIsaStatCols,
        clustermapsFormats: clustermapsFormats,
        clustermapsFigWidth: clustermapsFigWidth,
        clustermapsRowHeight: clustermapsRowHeight,
        clustermapsMinHeight: clustermapsMinHeight,
        clustermapsMaxHeight: clustermapsMaxHeight,
        clustermapsMitoSampleMode: clustermapsMitoSampleMode,
        clustermapsIsaAutoCandidates: clustermapsIsaAutoCandidates,
        spieceasiOutputDirAbs: spieceasiOutputDirAbs,
        spieceasiPrefix: spieceasiPrefix,
        spieceasiMinRelAbund: spieceasiMinRelAbund,
        spieceasiMinPrevalence: spieceasiMinPrevalence,
        spieceasiMethod: spieceasiMethod,
        spieceasiLambdaMinRatio: spieceasiLambdaMinRatio,
        spieceasiNlambda: spieceasiNlambda,
        spieceasiRepNum: spieceasiRepNum,
        spieceasiThresh: spieceasiThresh,
        spieceasiPulsarCriterion: spieceasiPulsarCriterion,
        spieceasiNcores: spieceasiNcores,
        spieceasiSeed: spieceasiSeed,
        spieceasiEdgeThreshold: spieceasiEdgeThreshold,
        spieceasiLayoutIters: spieceasiLayoutIters,
        networkGraphAllPath: networkGraphAllPath,
        networkGraphThrPath: networkGraphThrPath,
        networkNodeFeaturesPath: networkNodeFeaturesPath,
        networkLayoutSeed: networkLayoutSeed,
        networkLayoutScale: networkLayoutScale,
        networkDegreeScale: networkDegreeScale,
        networkDegreeSizeMode: networkDegreeSizeMode,
        networkDegreeMinArea: networkDegreeMinArea,
        networkEdgeWidthScale: networkEdgeWidthScale,
        networkIsaScale: networkIsaScale,
        networkAbundanceSizeMode: networkAbundanceSizeMode,
        networkAbundanceReference: networkAbundanceReference,
        networkAbundanceReferenceArea: networkAbundanceReferenceArea,
        networkAbundanceMinArea: networkAbundanceMinArea,
        networkAbundanceMaxArea: networkAbundanceMaxArea,
        networkAbundanceScalePower: networkAbundanceScalePower,
        networkMaxLabels: networkMaxLabels,
        networkModuleBestMinSize: networkModuleBestMinSize,
        networkModuleBestMinStability: networkModuleBestMinStability,
        networkModuleBestTopN: networkModuleBestTopN,
        networkModuleIsaSource: networkModuleIsaSource,
        networkModuleIsaMinStat: networkModuleIsaMinStat,
        networkModuleIsaMaxQ: networkModuleIsaMaxQ,
        networkMetadataPath: networkMetadataPath,
        networkColorCol: networkColorCol,
        networkIsaOverlayGroupsCsv: networkIsaOverlayGroupsCsv,
        networkGroupPaletteJson: networkGroupPaletteJson,
        networkGroupOrderJson: networkGroupOrderJson,
        networkFocusLabelJson: networkFocusLabelJson,
        networkModulePrimaryMethod: networkModulePrimaryMethod,
        networkModuleReps: networkModuleReps,
        networkModuleConsensusThreshold: networkModuleConsensusThreshold,
        networkModuleSeed: networkModuleSeed,
        networkModulesSubPath: networkModulesSubPath,
        networkModulesAllPath: networkModulesAllPath,
        masterSummaryOutputDirAbs: masterSummaryOutputDirAbs,
        masterSummaryClustermapsDirAbs: masterSummaryClustermapsDirAbs,
        masterSummaryIndicspeciesDirAbs: masterSummaryIndicspeciesDirAbs,
        masterSummarySpieceasiDirAbs: masterSummarySpieceasiDirAbs,
        masterSummaryAsvMagDirAbs: masterSummaryAsvMagDirAbs,
        masterSummaryWhitelistCsv: masterSummaryWhitelistCsv,
        masterSummaryMaxDirectCols: masterSummaryMaxDirectCols,
        asvMagLinkBarrnapDir: asvMagLinkBarrnapDir,
        asvMagLinkMasterTsv: asvMagLinkMasterTsv,
        asvMagLinkGenomeDir: asvMagLinkGenomeDir,
        asvMagLinkGenomeQcDir: asvMagLinkGenomeQcDir,
        asvMagLinkAutoBarrnap: asvMagLinkAutoBarrnap,
        asvMagLinkOutputDirAbs: asvMagLinkOutputDirAbs,
        asvMagLinkThreads: asvMagLinkThreads,
        asvMagLinkMinPident: asvMagLinkMinPident,
        asvMagLinkMinQcov: asvMagLinkMinQcov,
        asvMagLinkAsvTaxonomyMinConfidence: asvMagLinkAsvTaxonomyMinConfidence,
        asvMagLinkMinCompleteness: asvMagLinkMinCompleteness,
        asvMagLinkMaxContamination: asvMagLinkMaxContamination,
        asvMagLinkGuncAssessmentValue: asvMagLinkGuncAssessmentValue,
        asvMagLinkRequireSpeciesAssignment: asvMagLinkRequireSpeciesAssignment,
        asvMagLinkMinRrnaMarkerCount: asvMagLinkMinRrnaMarkerCount,
        asvMagLinkTopN: asvMagLinkTopN,
        asvMagLinkPlotTopN: asvMagLinkPlotTopN,
        asvMagNetworkOutputDirAbs: asvMagNetworkOutputDirAbs,
        asvMagNetworkPrefix: asvMagNetworkPrefix,
        asvMagNetworkAnchorTopN: asvMagNetworkAnchorTopN,
        asvMagNetworkGraphVariant: asvMagNetworkGraphVariant,
        asvMagNetworkMinPident: asvMagNetworkMinPident,
        asvMagNetworkMinQcov: asvMagNetworkMinQcov,
        asvMagNetworkAsvTaxonomySource: asvMagNetworkAsvTaxonomySource,
        asvMagNetworkMagTaxonomySource: asvMagNetworkMagTaxonomySource,
        asvMagNetworkMagAbundance: asvMagNetworkMagAbundance,
        asvMagNetworkMagIdMode: asvMagNetworkMagIdMode,
        asvMagNetworkMagAbundanceFormat: asvMagNetworkMagAbundanceFormat,
        asvMagNetworkMagAbundanceGenomeCol: asvMagNetworkMagAbundanceGenomeCol,
        asvMagNetworkMagAbundanceSampleCol: asvMagNetworkMagAbundanceSampleCol,
        asvMagNetworkMagAbundanceValueCol: asvMagNetworkMagAbundanceValueCol,
        asvMagNetworkMagAbundanceSeqkit: asvMagNetworkMagAbundanceSeqkit,
        asvMagNetworkMagAbundanceSeqkitFileCol: asvMagNetworkMagAbundanceSeqkitFileCol,
        asvMagNetworkMagAbundanceSeqkitCountCol: asvMagNetworkMagAbundanceSeqkitCountCol,
        asvMagNetworkMagTranscriptAbundance: asvMagNetworkMagTranscriptAbundance,
        asvMagNetworkMagTranscriptAbundanceFormat: asvMagNetworkMagTranscriptAbundanceFormat,
        asvMagNetworkMagTranscriptAbundanceGenomeCol: asvMagNetworkMagTranscriptAbundanceGenomeCol,
        asvMagNetworkMagTranscriptAbundanceSampleCol: asvMagNetworkMagTranscriptAbundanceSampleCol,
        asvMagNetworkMagTranscriptAbundanceValueCol: asvMagNetworkMagTranscriptAbundanceValueCol,
        asvMagNetworkMagTranscriptAbundanceSeqkit: asvMagNetworkMagTranscriptAbundanceSeqkit,
        asvMagNetworkMagTranscriptAbundanceSeqkitFileCol: asvMagNetworkMagTranscriptAbundanceSeqkitFileCol,
        asvMagNetworkMagTranscriptAbundanceSeqkitCountCol: asvMagNetworkMagTranscriptAbundanceSeqkitCountCol,
        asvMagNetworkMagTranscriptAbundanceNormalization: asvMagNetworkMagTranscriptAbundanceNormalization,
        asvMagCurtainsEnabled: asvMagCurtainsEnabled,
        asvMagCurtainsOutputDirAbs: asvMagCurtainsOutputDirAbs,
        asvMagCurtainsFormats: asvMagCurtainsFormats,
        asvMagCurtainsMaximumDepth: asvMagCurtainsMaximumDepth,
        asvMagCurtainsMinimumPointArea: asvMagCurtainsMinimumPointArea,
        asvMagCurtainsMaximumPointArea: asvMagCurtainsMaximumPointArea,
        asvMagCurtainsAsvSamplePointArea: asvMagCurtainsAsvSamplePointArea,
        asvMagCurtainsContrastPseudocountFpm: asvMagCurtainsContrastPseudocountFpm,
        asvMagCurtainsContrastFoldThreshold: asvMagCurtainsContrastFoldThreshold,
        asvMagCurtainsRenewalEvents: asvMagCurtainsRenewalEvents,
        asvMagCurtainsRenewalDateCol: asvMagCurtainsRenewalDateCol,
        asvMagNetworkMinSharedSamples: asvMagNetworkMinSharedSamples,
        asvMagNetworkMagKnn: asvMagNetworkMagKnn,
        asvMagNetworkAbundanceTransform: asvMagNetworkAbundanceTransform,
        asvMagNetworkMagAbundanceNormalization: asvMagNetworkMagAbundanceNormalization,
        asvMagNetworkAmbiguityTargetModules: asvMagNetworkAmbiguityTargetModules,
        asvMagNetworkAmbiguityFormats: asvMagNetworkAmbiguityFormats,
        asvMagNetworkFontFamily: asvMagNetworkFontFamily,
        asvMagNetworkBiochemicalGroupings: asvMagNetworkBiochemicalGroupings,
        asvMagNetworkIsaQThreshold: asvMagNetworkIsaQThreshold,
        asvMagNetworkGroupingSampleCol: asvMagNetworkGroupingSampleCol,
        asvMagNetworkGroupingCruiseCol: asvMagNetworkGroupingCruiseCol,
        asvMagNetworkGroupingGroups: asvMagNetworkGroupingGroups,
        asvMagNetworkGroupingCruiseGroups: asvMagNetworkGroupingCruiseGroups,
        asvMagNetworkGroupingPermutations: asvMagNetworkGroupingPermutations,
        groupGuildFunctionOutputDirAbs: groupGuildFunctionOutputDirAbs,
        groupGuildFunctionGroupings: groupGuildFunctionGroupings,
        groupGuildFunctionPrimary: groupGuildFunctionPrimary,
        groupGuildFunctionTargets: groupGuildFunctionTargets,
        groupGuildFunctionSampleCol: groupGuildFunctionSampleCol,
        groupGuildFunctionCruiseCol: groupGuildFunctionCruiseCol,
        groupGuildFunctionPermutations: groupGuildFunctionPermutations,
        groupGuildFunctionSeed: groupGuildFunctionSeed,
        groupGuildFunctionFormats: groupGuildFunctionFormats,
        groupGuildFunctionFontFamily: groupGuildFunctionFontFamily,
        groupGuildFunctionStrictFont: groupGuildFunctionStrictFont,
        groupGuildFunctionMaxModules: groupGuildFunctionMaxModules,
        groupGuildFunctionHeatmapTopAsvs: groupGuildFunctionHeatmapTopAsvs,
        groupGuildFunctionPcaLabelMinMeanAbundancePct: groupGuildFunctionPcaLabelMinMeanAbundancePct,
        groupGuildFunctionPcaOverlayEnabled: groupGuildFunctionPcaOverlayEnabled,
        groupGuildFunctionPcaScores: groupGuildFunctionPcaScores,
        groupGuildFunctionPcaExplained: groupGuildFunctionPcaExplained,
        groupGuildFunctionHybridAssignments: groupGuildFunctionHybridAssignments,
        groupGuildFunctionHybridCentroids: groupGuildFunctionHybridCentroids,
        asvMagNetworkFunctionalModuleMinFraction: asvMagNetworkFunctionalModuleMinFraction,
        asvMagNetworkFunctionalAnnotations: asvMagNetworkFunctionalAnnotations,
        powerAnalysisOutputDirAbs: powerAnalysisOutputDirAbs,
        powerAnalysisSampleCol: powerAnalysisSampleCol,
        powerAnalysisPatientCol: powerAnalysisPatientCol,
        powerAnalysisCaseCol: powerAnalysisCaseCol,
        powerAnalysisTypeCol: powerAnalysisTypeCol,
        powerAnalysisSampleSizesCancer: powerAnalysisSampleSizesCancer,
        powerAnalysisSampleSizesStype: powerAnalysisSampleSizesStype,
        powerAnalysisNSimulations: powerAnalysisNSimulations,
        powerAnalysisNPerm: powerAnalysisNPerm,
        powerAnalysisAlpha: powerAnalysisAlpha,
        powerAnalysisSeed: powerAnalysisSeed,
        powerAnalysisTransform: powerAnalysisTransform,
        powerAnalysisContralateralTypes: powerAnalysisContralateralTypes,
        powerAnalysisIndicspeciesDir: powerAnalysisIndicspeciesDir,
        taxonomyPatientAwareOutputDirAbs: taxonomyPatientAwareOutputDirAbs,
        taxonomyPatientAwareSampleCol: taxonomyPatientAwareSampleCol,
        taxonomyPatientAwarePatientCol: taxonomyPatientAwarePatientCol,
        taxonomyPatientAwareCaseCol: taxonomyPatientAwareCaseCol,
        taxonomyPatientAwareTypeCol: taxonomyPatientAwareTypeCol,
        taxonomyPatientAwareCountCol: taxonomyPatientAwareCountCol,
        taxonomyPatientAwareTaxLevels: taxonomyPatientAwareTaxLevels,
        taxonomyPatientAwareSampleTypes: taxonomyPatientAwareSampleTypes,
        taxonomyPatientAwareComparisonGroups: taxonomyPatientAwareComparisonGroups,
        taxonomyPatientAwareRunComparison: taxonomyPatientAwareRunComparison,
        taxonomyPatientAwareMinPrevalence: taxonomyPatientAwareMinPrevalence,
        taxonomyPatientAwareContralateralCol: taxonomyPatientAwareContralateralCol,
        taxonomyPatientAwareCancerSiteCol: taxonomyPatientAwareCancerSiteCol,
        taxonomyPatientAwareLungSideCol: taxonomyPatientAwareLungSideCol,
        taxonomyPatientAwareContralateralValue: taxonomyPatientAwareContralateralValue,
        taxonomyPatientAwareContralateralTypes: taxonomyPatientAwareContralateralTypes,
        taxonomyPatientAwareTransform: taxonomyPatientAwareTransform,
        taxonomyPatientAwareAlpha: taxonomyPatientAwareAlpha,
        taxonomyPatientAwareTopN: taxonomyPatientAwareTopN,
        taxonomyPatientAwareTypePalette: taxonomyPatientAwareTypePalette,
        taxonomyPatientAwareCasePalette: taxonomyPatientAwareCasePalette,
        lungStatusAnalysisOutputDirAbs: lungStatusAnalysisOutputDirAbs,
        lungStatusAnalysisSampleCol: lungStatusAnalysisSampleCol,
        lungStatusAnalysisTypeCol: lungStatusAnalysisTypeCol,
        lungStatusAnalysisSampleTypes: lungStatusAnalysisSampleTypes,
        lungStatusAnalysisCaseCol: lungStatusAnalysisCaseCol,
        lungStatusAnalysisPatientCol: lungStatusAnalysisPatientCol,
        lungStatusAnalysisCancerSiteCol: lungStatusAnalysisCancerSiteCol,
        lungStatusAnalysisLungCodeCol: lungStatusAnalysisLungCodeCol,
        lungStatusAnalysisTumorSideCol: lungStatusAnalysisTumorSideCol,
        lungStatusAnalysisContralateralCol: lungStatusAnalysisContralateralCol,
        lungStatusAnalysisHealthyCol: lungStatusAnalysisHealthyCol,
        lungStatusAnalysisStatusCol: lungStatusAnalysisStatusCol,
        lungStatusAnalysisStatusAValue: lungStatusAnalysisStatusAValue,
        lungStatusAnalysisStatusBValue: lungStatusAnalysisStatusBValue,
        lungStatusAnalysisReferenceStatusValue: lungStatusAnalysisReferenceStatusValue,
        lungStatusAnalysisPermutations: lungStatusAnalysisPermutations,
        lungStatusAnalysisSeed: lungStatusAnalysisSeed,
        indicspeciesGroupCols: indicspeciesGroupCols,
        indicspeciesEnabled: indicspeciesEnabled,
        indicspeciesPlotEnabled: indicspeciesPlotEnabled,
        indicspeciesLabelFocusedAsvs: indicspeciesLabelFocusedAsvs,
        indicspeciesAlignedEnabled: indicspeciesAlignedEnabled,
        clustermapsEnabled: clustermapsEnabled,
        clustermapsExternalIsaTable: clustermapsExternalIsaTable,
        clustermapsGroup1Order: clustermapsGroup1Order,
        clustermapsRunMito: clustermapsRunMito,
        spieceasiEnabled: spieceasiEnabled,
        spieceasiTranspose: spieceasiTranspose,
        spieceasiRemoveZeroVar: spieceasiRemoveZeroVar,
        spieceasiKeepNegative: spieceasiKeepNegative,
        spieceasiAllPosOnly: spieceasiAllPosOnly,
        spieceasiForceFilter: spieceasiForceFilter,
        spieceasiForceSpieceasi: spieceasiForceSpieceasi,
        spieceasiForceGraphs: spieceasiForceGraphs,
        spieceasiForceKeepIsaAsvs: spieceasiForceKeepIsaAsvs,
        spieceasiForceKeepAsvMagAsvs: spieceasiForceKeepAsvMagAsvs,
        networkEnabled: networkEnabled,
        networkModes: networkModes,
        networkModuleBestOnly: networkModuleBestOnly,
        networkModuleIsaOnly: networkModuleIsaOnly,
        networkModuleColorByIsa: networkModuleColorByIsa,
        networkModuleSubnetworksEnabled: networkModuleSubnetworksEnabled,
        networkModuleSubnetworkProminenceMetrics: networkModuleSubnetworkProminenceMetrics,
        networkModuleSubnetworkProminenceThreshold: networkModuleSubnetworkProminenceThreshold,
        networkModuleSubnetworkLabelTopN: networkModuleSubnetworkLabelTopN,
        networkModuleSubnetworkProminenceMinArea: networkModuleSubnetworkProminenceMinArea,
        networkModuleSubnetworkProminenceMaxArea: networkModuleSubnetworkProminenceMaxArea,
        networkModulesEnabled: networkModulesEnabled,
        genomeCooccurrenceEnabled: genomeCooccurrenceEnabled,
        genomeCooccurrenceOutputDirAbs: genomeCooccurrenceOutputDirAbs,
        genomeCooccurrenceGenomeQc: genomeCooccurrenceGenomeQc,
        genomeCooccurrenceMetagenome: genomeCooccurrenceMetagenome,
        genomeCooccurrenceMetatranscriptome: genomeCooccurrenceMetatranscriptome,
        genomeCooccurrenceMetadata: genomeCooccurrenceMetadata,
        genomeCooccurrenceSampleCol: genomeCooccurrenceSampleCol,
        genomeCooccurrenceGenomeCol: genomeCooccurrenceGenomeCol,
        genomeCooccurrenceMetagenomeValueCol: genomeCooccurrenceMetagenomeValueCol,
        genomeCooccurrenceMetatranscriptomeValueCol: genomeCooccurrenceMetatranscriptomeValueCol,
        genomeCooccurrenceInputScale: genomeCooccurrenceInputScale,
        genomeCooccurrenceInputFeatureLevel: genomeCooccurrenceInputFeatureLevel,
        genomeCooccurrenceExpectedSpecies: genomeCooccurrenceExpectedSpecies,
        genomeCooccurrenceClosureTotal: genomeCooccurrenceClosureTotal,
        genomeCooccurrenceMinimumReadCount: genomeCooccurrenceMinimumReadCount,
        genomeCooccurrenceMetadataSampleCol: genomeCooccurrenceMetadataSampleCol,
        genomeCooccurrenceCruiseCol: genomeCooccurrenceCruiseCol,
        genomeCooccurrenceSeasonCol: genomeCooccurrenceSeasonCol,
        genomeCooccurrenceDepthCol: genomeCooccurrenceDepthCol,
        genomeCooccurrenceMonthCol: genomeCooccurrenceMonthCol,
        genomeCooccurrenceSampleIdRegex: genomeCooccurrenceSampleIdRegex,
        genomeCooccurrenceMinRelAbund: genomeCooccurrenceMinRelAbund,
        genomeCooccurrenceMinPrevalence: genomeCooccurrenceMinPrevalence,
        genomeCooccurrenceZeroReplacementFraction: genomeCooccurrenceZeroReplacementFraction,
        genomeCooccurrenceMinAbsRho: genomeCooccurrenceMinAbsRho,
        genomeCooccurrenceBootstrapIterations: genomeCooccurrenceBootstrapIterations,
        genomeCooccurrenceMinBootstrapRecovery: genomeCooccurrenceMinBootstrapRecovery,
        genomeCooccurrenceMinSignConsistency: genomeCooccurrenceMinSignConsistency,
        genomeCooccurrencePermutations: genomeCooccurrencePermutations,
        genomeCooccurrencePermutationStrata: genomeCooccurrencePermutationStrata,
        genomeCooccurrenceMaxQ: genomeCooccurrenceMaxQ,
        genomeCooccurrenceSeed: genomeCooccurrenceSeed,
        networkModuleMethods: networkModuleMethods,
        networkModuleResolutions: networkModuleResolutions,
        masterSummaryEnabled: masterSummaryEnabled,
        asvMagLinkEnabled: asvMagLinkEnabled,
        asvMagNetworkEnabled: asvMagNetworkEnabled,
        groupGuildFunctionEnabled: groupGuildFunctionEnabled,
        asvMagLinkGenomeQcDirs: asvMagLinkGenomeQcDirs,
        asvMagLinkIdTokenIndexes: asvMagLinkIdTokenIndexes,
        powerAnalysisEnabled: powerAnalysisEnabled,
        powerAnalysisSkipEstimate: powerAnalysisSkipEstimate,
        powerAnalysisSkipPlot: powerAnalysisSkipPlot,
        powerAnalysisKeepContralateralInCancer: powerAnalysisKeepContralateralInCancer,
        taxonomyPatientAwareEnabled: taxonomyPatientAwareEnabled,
        taxonomyPatientAwareExcludeContralateral: taxonomyPatientAwareExcludeContralateral,
        taxonomyPatientAwareSkipOmnibus: taxonomyPatientAwareSkipOmnibus,
        lungStatusAnalysisEnabled: lungStatusAnalysisEnabled,
    ]
}

workflow ASPIRE_PREPROCESS {
    main:
    def rawReadsForAsv = rawReadsChannel
    def fastp_result = FASTP_QC(rawReadsForAsv)
    def reads_after_qc = fastp_result.reads
    def reads_after_merge = MERGE_READS(reads_after_qc)
    def reads_after_filter = FILTER_READS(reads_after_merge)
    def reads_for_concat = reads_after_filter
    if( concatRelabelEnabled ) {
        def relabeled_stage = RELABEL_FILTERED(reads_after_filter)
        reads_for_concat = relabeled_stage.relabeled
    }
    def relabeled_fasta_files = reads_for_concat.map { parts -> parts[1] }
    def concat_stage = CONCAT_FASTAS(relabeled_fasta_files.collect())
    def concat_for_derep = concat_stage.concat_for_derep
    def concat_for_counts = concat_stage.concat_for_counts

    def derep_input = DEREPLICATE(concat_for_derep)
    def denoise_input = DENOISE(derep_input)
    def nochi_input = CHIMERA_CHECK(denoise_input)

    def count_matrix_stage = CREATE_COUNT_MATRIX(concat_for_counts, nochi_input)
    def count_matrix_channel = count_matrix_stage.count_matrix
    def asv_counts_for_sankey = count_matrix_channel.map { tuple -> tuple[0] }
    def filtered_stage = FILTER_TABLE(count_matrix_channel)
    def filtered_channel = filtered_stage.filtered
    def filtered_fasta_for_taxonomy = filtered_channel.map { tuple -> tuple[1] }
    def sina_stage = SINA_TRIM(filtered_fasta_for_taxonomy)
    def taxonomy_stage = TAXONOMY(sina_stage.trimmed_fasta)
    def runMitoStages = mitoEnabled || filterCountsEnabled
    def filter_counts_stage = null
    if( runMitoStages ) {
        def blast_database_stage = PREPARE_BLAST_DATABASES(
            mitoBlastFastaPath ?: mitoBlastDbPath,
            mitoBiofFastaPath ?: mitoBiofDbPath
        )
        def mitomaster_stage = MITOMASTER(filtered_channel, blast_database_stage.databases)
        def mito_summary = MITO_DECONTAM(mitomaster_stage.mito_artifacts, taxonomy_stage.taxonomy_table)
        if( filterCountsEnabled ) {
            filter_counts_stage = FILTER_COUNTS(filtered_channel, taxonomy_stage.taxonomy_table, mito_summary.nontarget_table)
        }
    }
    if( metadataPlotsEnabled && filter_counts_stage == null ) {
        exit 1, "metadata_plots.enabled requires filter_counts outputs but filter_counts stage was not executed"
    }
    def general_stats_stage = null
    if( generalStatsEnabled ) {
        general_stats_stage = GENERAL_STATS(concat_for_counts)
    }

    if( filter_counts_stage == null || general_stats_stage == null ) {
        exit 1, "Canonical preprocessing export requires filter_counts.enabled=true and general_stats.enabled=true"
    }
    def preprocessing_dataset_stage = PREPROCESS_DATASET(
        count_matrix_channel.map { tuple -> tuple[0] },
        filtered_channel.map { tuple -> tuple[0] },
        filtered_channel.map { tuple -> tuple[1] },
        filter_counts_stage.filtered_counts,
        filter_counts_stage.filtered_micro,
        filter_counts_stage.filtered_mito,
        filter_counts_stage.filtered_decon,
        taxonomy_stage.taxonomy_table,
        general_stats_stage.fastq_stats,
        general_stats_stage.fastp_stats,
        general_stats_stage.filtered_stats,
        general_stats_stage.concat_stats,
        file(manifestPath, checkIfExists: true),
        preprocessingParametersJson
    )

    emit:
    dataset = preprocessing_dataset_stage.dataset
}

workflow ASPIRE_ANALYSIS {
    take:
    preprocessing_dataset

    main:
    def validated_preprocessing = VALIDATE_PREPROCESS_DATASET(preprocessing_dataset)
    def datasetFile = { String name -> validated_preprocessing.dataset.map { directory -> file("${directory}/${name}") } }
    def count_matrix_channel = datasetFile('asv_counts.tsv')
    def asv_counts_for_sankey = count_matrix_channel.map { it }
    def filtered_fasta_for_taxonomy = datasetFile('asv_sequences.filtered.fasta.gz')
    def taxonomy_stage = [taxonomy_table: datasetFile('taxonomy.tsv')]
    def filter_counts_stage = [
        filtered_counts: datasetFile('asv_target.tsv'),
        filtered_micro: datasetFile('asv_target.micro.tsv'),
        filtered_mito: datasetFile('asv_target.mito.tsv'),
        filtered_decon: datasetFile('asv_target.decon.tsv')
    ]
    def general_stats_stage = [
        fastq_stats: datasetFile('fastq_stats.tsv'),
        fastp_stats: datasetFile('fastp_fastqs.tsv'),
        filtered_stats: datasetFile('filtered_fastas.tsv'),
        concat_stats: datasetFile('concat_fastas.tsv')
    ]
    def sample_manifest = datasetFile('sample_manifest.tsv')

    def metadata_analysis_stage = null
    def metaMicroForNetwork = null
    def asvFinalForSpieceasi = null
    def asvFinalForNetwork = null
    def asvMetaForMasterSummary = null
    def asvFinalForMasterSummary = null
    def indicspeciesTablesForOverlay = Channel.value(file(emptyModulesPath))
    def indicspeciesGroup1SummaryForSpieceasi = Channel.value(file(emptyModulesPath))
    if( metadataPlotsEnabled ) {
        metadata_analysis_stage = RUN_METADATA_ANALYSES(
            general_stats_stage.fastq_stats,
            filter_counts_stage.filtered_counts,
            filter_counts_stage.filtered_mito,
            taxonomy_stage.taxonomy_table,
            sample_manifest
        )
        metaMicroForNetwork = metadata_analysis_stage.meta_micro_network
        asvFinalForSpieceasi = metadata_analysis_stage.asv_final_spieceasi
        asvFinalForNetwork = metadata_analysis_stage.asv_final_network
        asvMetaForMasterSummary = metadata_analysis_stage.asv_meta_master_summary
        asvFinalForMasterSummary = metadata_analysis_stage.asv_final_master_summary
        if( indicspeciesRunEarly ) {
            indicspeciesTablesForOverlay = metadata_analysis_stage.indicspecies_tables
        }
        if( spieceasiForceKeepIsaAsvs ) {
            indicspeciesGroup1SummaryForSpieceasi = metadata_analysis_stage.indicspecies_group1_summary
        }
    }
    def asv_mag_link_stage = null
    if( asvMagLinkEnabled ) {
        asv_mag_link_stage = ASV_MAG_LINK(
            filtered_fasta_for_taxonomy,
            taxonomy_stage.taxonomy_table
        )
    }
    def asvMagLinksForSpieceasi = asv_mag_link_stage != null ?
        asv_mag_link_stage.pairing :
        Channel.value(file(emptyModulesPath))
    def spieceasi_stage = null
    def graphAllForModules = null
    def graphAllForNetwork = null
    def graphAllForAsvMagNetwork = null
    def graphAllForModuleMeasurement = null
    def graphThrForModules = null
    def graphThrForNetwork = null
    def graphThrForAsvMagNetwork = null
    def nodeFeaturesForNetwork = null
    def nodeFeaturesForAsvMagNetwork = null
    def nodeFeaturesForModuleMeasurement = null
    def modulesSubForNetwork = null
    def modulesAllForNetwork = null
    if( spieceasiEnabled ) {
        spieceasi_stage = SPIECEASI(
            asvFinalForSpieceasi,
            indicspeciesGroup1SummaryForSpieceasi,
            asvMagLinksForSpieceasi
        )
        graphAllForModules = spieceasi_stage.graph_all.map { it }
        graphAllForNetwork = spieceasi_stage.graph_all.map { it }
        graphAllForAsvMagNetwork = spieceasi_stage.graph_all.map { it }
        graphAllForModuleMeasurement = spieceasi_stage.graph_all.map { it }
        graphThrForModules = spieceasiAllPosOnly ? spieceasi_stage.graph_all.map { it } : spieceasi_stage.graph_thr.map { it }
        graphThrForNetwork = spieceasiAllPosOnly ? spieceasi_stage.graph_all.map { it } : spieceasi_stage.graph_thr.map { it }
        graphThrForAsvMagNetwork = spieceasiAllPosOnly ? spieceasi_stage.graph_all.map { it } : spieceasi_stage.graph_thr.map { it }
        nodeFeaturesForNetwork = spieceasi_stage.node_features
        nodeFeaturesForAsvMagNetwork = spieceasi_stage.node_features.map { it }
        nodeFeaturesForModuleMeasurement = spieceasi_stage.node_features.map { it }
    } else if( networkEnabled ) {
        graphAllForModules = Channel.value(file(networkGraphAllPath))
        graphAllForNetwork = Channel.value(file(networkGraphAllPath))
        graphAllForAsvMagNetwork = Channel.value(file(networkGraphAllPath))
        graphAllForModuleMeasurement = Channel.value(file(networkGraphAllPath))
        graphThrForModules = Channel.value(file(spieceasiAllPosOnly ? networkGraphAllPath : networkGraphThrPath))
        graphThrForNetwork = Channel.value(file(spieceasiAllPosOnly ? networkGraphAllPath : networkGraphThrPath))
        graphThrForAsvMagNetwork = Channel.value(file(spieceasiAllPosOnly ? networkGraphAllPath : networkGraphThrPath))
        nodeFeaturesForNetwork = Channel.value(file(networkNodeFeaturesPath))
        nodeFeaturesForAsvMagNetwork = Channel.value(file(networkNodeFeaturesPath))
        nodeFeaturesForModuleMeasurement = Channel.value(file(networkNodeFeaturesPath))
    }
    def measurement_association_stage = null
    if( measurementAssociationEnabled ) {
        def measurementAssociationSubsetAudit = Channel.value(file(emptyModulesPath))
        if( measurementAssociationSubsetSource == 'spieceasi_standard_filter' ) {
            if( spieceasi_stage == null ) {
                exit 1, "measurement_association.asv_subset_source=spieceasi_standard_filter requires spieceasi.enabled=true"
            }
            measurementAssociationSubsetAudit = spieceasi_stage.filter_audit.map { it }
        }
        measurement_association_stage = MEASUREMENT_ASSOCIATION(
            metadata_analysis_stage.asv_meta_measurement_association,
            metadata_analysis_stage.meta_micro_measurement_association,
            metadata_analysis_stage.asv_final_measurement_association,
            measurementAssociationSubsetAudit
        )
    }
    def titan_collect_stage = null
    if( titanEnabled ) {
        if( spieceasi_stage == null || measurement_association_stage == null ) {
            exit 1, "TITAN requires completed SPIEC-EASI and measurement-association stages"
        }
        def titan_prepare_stage = TITAN_PREPARE(
            asvFinalForSpieceasi.map { it },
            spieceasi_stage.filter_audit.map { it },
            taxonomy_stage.taxonomy_table.map { it },
            measurement_association_stage.done.map { it }
        )
        def titan_install_stage = TITAN_INSTALL()
        def titan_analysis_stage = TITAN_ANALYSIS(
            titan_prepare_stage.prepared,
            titan_install_stage.library
        )
        titan_collect_stage = TITAN_COLLECT(
            titan_prepare_stage.prepared.map { it },
            titan_analysis_stage.raw,
            titan_install_stage.manifest
        )
        if( titanPlotEnabled ) {
            TITAN_PLOTS(titan_collect_stage.done)
        }
    }
    def microbial_compartment_posthoc_stage = null
    if( microbialCompartmentEnabled ) {
        if( spieceasi_stage == null || metadata_analysis_stage == null ) {
            exit 1, "Microbial-compartment inference requires completed SPIEC-EASI and metadata stages"
        }
        def microbial_compartment_prepare_stage = MICROBIAL_COMPARTMENT_PREPARE(
            asvFinalForSpieceasi.map { it },
            spieceasi_stage.filter_audit.map { it },
            metadata_analysis_stage.meta_micro_network.map { it }
        )
        def microbial_compartment_inference_stage = MICROBIAL_COMPARTMENT_INFERENCE(
            microbial_compartment_prepare_stage.prepared
        )
        microbial_compartment_posthoc_stage = MICROBIAL_COMPARTMENT_POSTHOC(
            microbial_compartment_prepare_stage.prepared.map { it },
            microbial_compartment_inference_stage.inference,
            metadata_analysis_stage.meta_micro_network.map { it },
            taxonomy_stage.taxonomy_table.map { it }
        )
        if( microbialCompartmentPlotEnabled ) {
            MICROBIAL_COMPARTMENT_PLOTS(microbial_compartment_posthoc_stage.done)
        }
        if( indicspeciesRequiresMicrobialCompartments ) {
            def microbial_compartment_indicspecies_stage = INDICSPECIES(
                microbial_compartment_posthoc_stage.metadata,
                asvFinalForSpieceasi.map { it }
            )
            indicspeciesTablesForOverlay = microbial_compartment_indicspecies_stage.all_tables.collect()
            if( indicspeciesPlotEnabled ) {
                INDICSPECIES_PLOTS(
                    microbial_compartment_posthoc_stage.metadata.map { it },
                    microbial_compartment_indicspecies_stage.all_tables.collect()
                )
            }
            if( indicspeciesAlignedEnabled ) {
                INDICSPECIES_ALIGNED_PLOTS(
                    microbial_compartment_indicspecies_stage.all_tables.collect()
                )
            }
        }
        if( groupingDiagnosticsRequiresMicrobialCompartments ) {
            GROUPING_DIAGNOSTICS(
                microbial_compartment_posthoc_stage.metadata.map { it },
                asvFinalForSpieceasi.map { it }
            )
        }
    }
    def community_turnover_lcbd_stage = null
    if( communityTurnoverEnabled ) {
        if( spieceasi_stage == null || microbial_compartment_posthoc_stage == null ) {
            exit 1, "Community turnover requires completed SPIEC-EASI and microbial-compartment metadata stages"
        }
        def community_turnover_prepare_stage = COMMUNITY_TURNOVER_PREPARE(
            asvFinalForSpieceasi.map { it },
            spieceasi_stage.filter_audit.map { it },
            microbial_compartment_posthoc_stage.metadata.map { it }
        )
        def community_turnover_analysis_stage = COMMUNITY_TURNOVER_ANALYSIS(
            community_turnover_prepare_stage.prepared,
            microbial_compartment_posthoc_stage.metadata.map { it }
        )
        community_turnover_lcbd_stage = COMMUNITY_TURNOVER_LCBD(
            community_turnover_prepare_stage.prepared.map { it },
            microbial_compartment_posthoc_stage.metadata.map { it },
            community_turnover_analysis_stage.analysis
        )
        if( communityTurnoverPlotEnabled ) {
            COMMUNITY_TURNOVER_PLOTS(community_turnover_lcbd_stage.done)
        }
    }
    if( networkEnabled ) {
        if( networkModulesEnabled ) {
            def network_modules_stage = NETWORK_MODULES(
                graphAllForModules,
                graphThrForModules
            )
            modulesSubForNetwork = network_modules_stage.modules_sub
            modulesAllForNetwork = network_modules_stage.modules_all
        } else {
            def modulesSubFile = new File(networkModulesSubPath)
            def modulesAllFile = new File(networkModulesAllPath)
            modulesSubForNetwork = modulesSubFile.exists() ? Channel.value(file(networkModulesSubPath)) : Channel.value(file(emptyModulesPath))
            modulesAllForNetwork = modulesAllFile.exists() ? Channel.value(file(networkModulesAllPath)) : Channel.value(file(emptyModulesPath))
        }
    }
    def microbial_state_interpretation_stage = null
    if( microbialStateInterpretationEnabled ) {
        if( microbial_compartment_posthoc_stage == null || modulesAllForNetwork == null || measurement_association_stage == null ) {
            exit 1, "Microbial-state interpretation requires completed microbial-compartment, network-module, and measurement-association stages"
        }
        microbial_state_interpretation_stage = MICROBIAL_STATE_INTERPRETATION(
            microbial_compartment_posthoc_stage.done,
            modulesAllForNetwork.map { it },
            microbial_compartment_posthoc_stage.metadata.map { it },
            taxonomy_stage.taxonomy_table.map { it },
            measurement_association_stage.done,
            ((config.microbial_state_interpretation?.contour_time_step_days ?: 7) as int),
            ((config.microbial_state_interpretation?.contour_max_time_support_days ?: 90.0d) as double),
            ((config.microbial_state_interpretation?.contour_max_depth_support_m ?: 30.0d) as double),
            (config.microbial_state_interpretation?.contour_reference_grid ?
                resolveOptionalPath(config.microbial_state_interpretation.contour_reference_grid.toString(), configRoot) : '')
        )
    }
    def graph_network_stage = null
    if( networkEnabled ) {
        def networkMetadataChannel = microbial_compartment_posthoc_stage != null ?
            microbial_compartment_posthoc_stage.metadata.map { it } :
            (metaMicroForNetwork != null ? metaMicroForNetwork : Channel.value(file(networkMetadataPath)))
        def asvMagReadyForNetwork = asv_mag_link_stage != null ? asv_mag_link_stage.done : Channel.value(file(emptyModulesPath))
        def microbialStateReadyForNetwork = microbial_state_interpretation_stage != null ?
            microbial_state_interpretation_stage.done : Channel.value(file(emptyModulesPath))
        graph_network_stage = RUN_GRAPH_NETWORK(
            graphAllForNetwork,
            graphThrForNetwork,
            nodeFeaturesForNetwork,
            asvFinalForNetwork,
            networkMetadataChannel,
            asvMagReadyForNetwork,
            taxonomy_stage.taxonomy_table,
            indicspeciesTablesForOverlay,
            modulesSubForNetwork,
            modulesAllForNetwork,
            microbialStateReadyForNetwork
        )
    }
    def module_measurement_association_stage = null
    if( moduleMeasurementAssociationEnabled ) {
        if( modulesAllForNetwork == null || measurement_association_stage == null || graph_network_stage == null ) {
            exit 1, "Module-measurement association requires network modules, network analysis, and the measurement-association stage"
        }
        module_measurement_association_stage = MODULE_MEASUREMENT_ASSOCIATION(
            metadata_analysis_stage.asv_final_measurement_association,
            filter_counts_stage.filtered_micro.map { it },
            microbial_compartment_posthoc_stage.metadata.map { it },
            modulesAllForNetwork.map { it },
            graphAllForModuleMeasurement,
            graph_network_stage.layout_all,
            nodeFeaturesForModuleMeasurement,
            taxonomy_stage.taxonomy_table,
            measurement_association_stage.done,
            moduleMeasurementNetworkMetricTopN
        )
    }
    def asv_mag_network_stage = null
    if( asvMagNetworkEnabled ) {
        def asvMagReadyForAsvMagNetwork = asv_mag_link_stage != null ? asv_mag_link_stage.done : Channel.value(file(emptyModulesPath))
        def graphForAsvMagNetwork = asvMagNetworkGraphVariant == 'all' ? graphAllForAsvMagNetwork : graphThrForAsvMagNetwork
        def cruiseMetadataForAsvMagNetwork = asvMagNetworkCruiseMetadata ?
            Channel.value(file(asvMagNetworkCruiseMetadata)) :
            Channel.value(file(emptyModulesPath))
        asv_mag_network_stage = RUN_ASV_MAG_NETWORK(
            graphForAsvMagNetwork,
            nodeFeaturesForAsvMagNetwork,
            modulesAllForNetwork,
            graph_network_stage.module_best_stats_all,
            taxonomy_stage.taxonomy_table,
            metadata_analysis_stage.asv_final_spieceasi,
            cruiseMetadataForAsvMagNetwork,
            microbial_compartment_posthoc_stage.metadata.map { it },
            indicspeciesTablesForOverlay,
            asvMagReadyForAsvMagNetwork
        )
    }
    def genome_cooccurrence_stage = null
    if( genomeCooccurrenceEnabled ) {
        if( asv_mag_network_stage == null || titan_collect_stage == null || modulesAllForNetwork == null || microbial_state_interpretation_stage == null ) {
            exit 1, "Genome co-occurrence requires completed ASV-MAG mappings, TITAN results, ASV ecological modules, and microbial-state interpretation"
        }
        genome_cooccurrence_stage = GENOME_COOCCURRENCE(
            Channel.value(file(genomeCooccurrenceGenomeQc)),
            Channel.value(file(genomeCooccurrenceMetagenome)),
            Channel.value(file(genomeCooccurrenceMetatranscriptome)),
            Channel.value(file(genomeCooccurrenceMetadata)),
            asv_mag_network_stage.analysis_mappings,
            modulesAllForNetwork.map { it },
            indicspeciesTablesForOverlay,
            titan_collect_stage.taxon_results,
            microbial_state_interpretation_stage.done
        )
    }
    def asv_mag_curtains_stage = null
    if( asvMagCurtainsEnabled ) {
        if( asv_mag_network_stage == null || microbial_state_interpretation_stage == null ) {
            exit 1, "ASV-MAG curtains require completed ASV-MAG network and microbial-state interpretation stages"
        }
        asv_mag_curtains_stage = ASV_MAG_CURTAINS(
            asv_mag_network_stage.done,
            microbial_state_interpretation_stage.done
        )
    }
    def module_mag_anchors_stage = null
    if( networkEnabled && asvMagLinkEnabled ) {
        def moduleMagGraphDone = graph_network_stage != null ? graph_network_stage.done : Channel.value(file(emptyModulesPath))
        def moduleMagAsvDone = asv_mag_link_stage != null ? asv_mag_link_stage.done : Channel.value(file(emptyModulesPath))
        module_mag_anchors_stage = RUN_MODULE_MAG_ANCHORS(
            modulesAllForNetwork,
            nodeFeaturesForNetwork,
            taxonomy_stage.taxonomy_table,
            metadata_analysis_stage.asv_final_spieceasi,
            microbial_compartment_posthoc_stage.metadata.map { it },
            asv_mag_network_stage.analysis_mappings,
            moduleMagAsvDone,
            moduleMagGraphDone
        )
    }
    def group_guild_function_stage = null
    if( groupGuildFunctionEnabled ) {
        def groupGuildPcaScoresInput = groupGuildFunctionPcaOverlayEnabled ?
            Channel.value(file(groupGuildFunctionPcaScores)) :
            Channel.value(file(emptyModulesPath))
        def groupGuildPcaExplainedInput = groupGuildFunctionPcaOverlayEnabled ?
            Channel.value(file(groupGuildFunctionPcaExplained)) :
            Channel.value(file(emptyModulesPath))
        def groupGuildHybridAssignmentsInput = groupGuildFunctionPcaOverlayEnabled ?
            Channel.value(file(groupGuildFunctionHybridAssignments)) :
            Channel.value(file(emptyModulesPath))
        def groupGuildHybridCentroidsInput = groupGuildFunctionPcaOverlayEnabled ?
            Channel.value(file(groupGuildFunctionHybridCentroids)) :
            Channel.value(file(emptyModulesPath))
        group_guild_function_stage = RUN_GROUP_GUILD_FUNCTION(
            microbial_compartment_posthoc_stage.metadata.map { it },
            metadata_analysis_stage.asv_final_spieceasi,
            modulesAllForNetwork,
            nodeFeaturesForNetwork,
            indicspeciesTablesForOverlay,
            asv_mag_network_stage.analysis_mappings,
            asv_mag_network_stage.functional_modules,
            groupGuildPcaScoresInput,
            groupGuildPcaExplainedInput,
            groupGuildHybridAssignmentsInput,
            groupGuildHybridCentroidsInput,
            taxonomy_stage.taxonomy_table
        )
    }
    def ecological_context_atlas_stage = null
    if( ecologicalContextAtlasEnabled ) {
        if( titan_collect_stage == null || microbial_compartment_posthoc_stage == null || module_measurement_association_stage == null || asv_mag_network_stage == null || graph_network_stage == null || group_guild_function_stage == null ) {
            exit 1, "Ecological-context atlas requires completed TITAN, microbial compartments, module-measurement, module-compartment, network, and ASV-MAG network stages"
        }
        ecological_context_atlas_stage = ECOLOGICAL_CONTEXT_ATLAS(
            microbial_compartment_posthoc_stage.metadata.map { it },
            microbial_compartment_posthoc_stage.done,
            modulesAllForNetwork.map { it },
            nodeFeaturesForNetwork.map { it },
            graph_network_stage.layout_all,
            graphAllForNetwork.map { it },
            taxonomy_stage.taxonomy_table.map { it },
            titan_collect_stage.done,
            measurement_association_stage.done,
            module_measurement_association_stage.done,
            group_guild_function_stage.done,
            asv_mag_network_stage.analysis_mappings,
            asv_mag_network_stage.functional_modules,
            asv_mag_network_stage.done
        )
    }
    def sankey_stage = null
    if( sankeyEnabled ) {
        sankey_stage = SANKEY(
            general_stats_stage.fastq_stats,
            general_stats_stage.filtered_stats,
            asv_counts_for_sankey,
            filter_counts_stage.filtered_decon,
            filter_counts_stage.filtered_micro,
            sample_manifest
        )
    }
    if( masterSummaryEnabled ) {
        if( asvMetaForMasterSummary == null || asvFinalForMasterSummary == null ) {
            exit 1, "master_summary.enabled requires ASV_meta and ASV_final inputs from metadata/batch stages"
        }
        def masterSummaryNetworkDone = graph_network_stage != null ? graph_network_stage.done : Channel.value(file(emptyModulesPath))
        def masterSummarySankeyDone = sankey_stage != null ? sankey_stage.done : Channel.value(file(emptyModulesPath))
        def masterSummaryAsvMagDone = asv_mag_link_stage != null ? asv_mag_link_stage.done : Channel.value(file(emptyModulesPath))
        def masterSummaryOptionalDone = genome_cooccurrence_stage != null ? genome_cooccurrence_stage.done : Channel.value(file(emptyModulesPath))
        MASTER_SUMMARY(
            asvMetaForMasterSummary,
            asvFinalForMasterSummary,
            masterSummaryNetworkDone,
            masterSummarySankeyDone,
            masterSummaryAsvMagDone,
            masterSummaryOptionalDone
        )
    }
}

workflow {
    if( aspirePhase == 'preprocess' ) {
        ASPIRE_PREPROCESS()
    } else if( aspirePhase == 'analysis' ) {
        def canonicalDatasetPath = new File(publicOutputDir, 'preprocessing_dataset').canonicalPath
        if( !new File(canonicalDatasetPath, 'manifest.json').isFile() ) {
            exit 1, "Canonical preprocessing dataset not found at ${canonicalDatasetPath}. Run --phase preprocess first."
        }
        def canonicalDataset = file(canonicalDatasetPath, checkIfExists: true)
        ASPIRE_ANALYSIS(Channel.value(canonicalDataset))
    } else {
        def preprocessing = ASPIRE_PREPROCESS()
        ASPIRE_ANALYSIS(preprocessing.dataset)
    }
}

workflow RUN_METADATA_ANALYSES {
    take:
    fastq_stats
    filtered_micro
    filtered_mito
    taxonomy_table
    sample_manifest

    main:
    metadata_stage = PLOT_METADATA(
        fastq_stats,
        filtered_micro,
        filtered_mito,
        taxonomy_table,
        sample_manifest
    )
    def baseMetadata = metadata_stage.metadata_micro
    def baseAsvMeta = metadata_stage.asv_meta_micro
    def baseAsvFinal = metadata_stage.asv_final_micro

    if( asvTimeDepthCurtainEnabled ) {
        ASV_TIME_DEPTH_CURTAIN(baseMetadata.map { it })
    }

    metaMicroForBatch = baseMetadata.map { it }
    metaMicroForOutlier = baseMetadata.map { it }
    metaMicroForPlotUpset = baseMetadata.map { it }
    metaMicroForCollectors = baseMetadata.map { it }
    metaMicroForDiversity = baseMetadata.map { it }
    metaMicroForIndicspecies = baseMetadata.map { it }
    metaMicroForIndicspeciesPlots = baseMetadata.map { it }
    metaMicroForClustermaps = baseMetadata.map { it }
    metaMicroForNetwork = baseMetadata.map { it }
    metaMicroForMeasurementAssociation = baseMetadata.map { it }
    metaMicroForGroupingDiagnostics = baseMetadata.map { it }
    metaMicroForCommunityPredictor = baseMetadata.map { it }
    asvMetaForBatch = baseAsvMeta.map { it }
    asvMetaForGroupAugmentation = baseAsvMeta.map { it }
    asvMetaSeedForCorrection = baseAsvMeta.map { it }
    asvMetaForBubbleplotter = baseAsvMeta.map { it }
    asvMetaForUmap = baseAsvMeta.map { it }
    asvMetaForClustermaps = baseAsvMeta.map { it }
    asvMetaForVocCorrelation = baseAsvMeta.map { it }
    asvMetaForMeasurementAssociation = baseAsvMeta.map { it }
    asvMetaForPowerAnalysis = baseAsvMeta.map { it }
    asvMetaForTaxonomyPatientAware = baseAsvMeta.map { it }
    asvMetaForLungStatus = baseAsvMeta.map { it }
    asvMetaForMasterSummary = baseAsvMeta.map { it }

    asvFinalForBatch = baseAsvFinal.map { it }
    asvFinalForCollectors = baseAsvFinal.map { it }
    asvFinalForDiversity = baseAsvFinal.map { it }
    asvFinalForIndicspecies = baseAsvFinal.map { it }
    asvFinalForSpieceasi = baseAsvFinal.map { it }
    asvFinalForNetwork = baseAsvFinal.map { it }
    asvFinalForVocCorrelation = baseAsvFinal.map { it }
    asvFinalForMeasurementAssociation = baseAsvFinal.map { it }
    asvFinalForGroupingDiagnostics = baseAsvFinal.map { it }
    asvFinalForCommunityPredictor = baseAsvFinal.map { it }
    asvFinalForPowerAnalysis = baseAsvFinal.map { it }
    asvFinalForTaxonomyPatientAware = baseAsvFinal.map { it }
    asvFinalForLungStatus = baseAsvFinal.map { it }
    asvFinalForMasterSummary = baseAsvFinal.map { it }

    grouping_diagnostics_stage = null
    if( groupingDiagnosticsRunEarly ) {
        grouping_diagnostics_stage = GROUPING_DIAGNOSTICS(
            metaMicroForGroupingDiagnostics,
            asvFinalForGroupingDiagnostics
        )
    }

    if( groupingDiagnosticsApplySoftLabels && groupingDiagnosticsRunEarly ) {
        group_label_augmentation_stage = GROUP_LABEL_AUGMENTATION(
            baseMetadata.map { it },
            asvMetaForGroupAugmentation,
            grouping_diagnostics_stage.soft_assignments,
            grouping_diagnostics_stage.soft_validation_summary
        )
        metaMicroForBatch = group_label_augmentation_stage.metadata_augmented.map { it }
        metaMicroForOutlier = group_label_augmentation_stage.metadata_augmented.map { it }
        metaMicroForPlotUpset = group_label_augmentation_stage.metadata_augmented.map { it }
        metaMicroForCollectors = group_label_augmentation_stage.metadata_augmented.map { it }
        metaMicroForDiversity = group_label_augmentation_stage.metadata_augmented.map { it }
        // Indicator analyses are publication-primary on observed BASIN labels.
        // Soft-label performance remains reported by GROUPING_DIAGNOSTICS as a
        // sensitivity analysis, but inferred labels must not define indicators.
        metaMicroForClustermaps = group_label_augmentation_stage.metadata_augmented.map { it }
        // Network ISA overlays and guild evidence inherit observed-label ISA.
        metaMicroForMeasurementAssociation = group_label_augmentation_stage.metadata_augmented.map { it }
        // Community predictor performs its own observed-label filtering.
        asvMetaForBatch = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaSeedForCorrection = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaForBubbleplotter = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaForUmap = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaForClustermaps = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaForVocCorrelation = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaForMeasurementAssociation = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaForPowerAnalysis = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaForTaxonomyPatientAware = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaForLungStatus = group_label_augmentation_stage.asv_meta_augmented.map { it }
        asvMetaForMasterSummary = group_label_augmentation_stage.asv_meta_augmented.map { it }
    }

    if( plotUpsetEnabled ) {
        PLOT_UPSET(metaMicroForPlotUpset)
    }

    batch_stage = null
    asvClrForOutlier = null
    if( batchCorrectionEnabled ) {
        batch_stage = ASV_BATCH_CORRECTION(
            metaMicroForBatch,
            asvMetaForBatch,
            asvFinalForBatch
        )
        asvClrForOutlier = batch_stage.asv_clr_selected
        asvFinalForCollectors = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForDiversity = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForIndicspecies = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForSpieceasi = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForNetwork = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForVocCorrelation = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForMeasurementAssociation = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForGroupingDiagnostics = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForCommunityPredictor = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForPowerAnalysis = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForTaxonomyPatientAware = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForLungStatus = batch_stage.asv_selected_counts_int.map { it }
        asvFinalForMasterSummary = batch_stage.asv_selected_counts_int.map { it }
        umapResultsForTrajectory = batch_stage.umap_results
        if( bubbleplotterEnabled || umapClusteringEnabled || clustermapsEnabled || vocCorrelationEnabled || measurementAssociationEnabled || masterSummaryEnabled || powerAnalysisEnabled || taxonomyPatientAwareEnabled || lungStatusAnalysisEnabled ) {
            corrected_asv_meta_stage = ASV_META_FROM_CORRECTED(
                asvMetaSeedForCorrection,
                batch_stage.asv_selected_counts_int
            )
            asvMetaForBubbleplotter = corrected_asv_meta_stage.asv_meta_corrected.map { it }
            asvMetaForUmap = corrected_asv_meta_stage.asv_meta_corrected.map { it }
            asvMetaForClustermaps = corrected_asv_meta_stage.asv_meta_corrected.map { it }
            asvMetaForVocCorrelation = corrected_asv_meta_stage.asv_meta_corrected.map { it }
            asvMetaForMeasurementAssociation = corrected_asv_meta_stage.asv_meta_corrected.map { it }
            asvMetaForPowerAnalysis = corrected_asv_meta_stage.asv_meta_corrected.map { it }
            asvMetaForTaxonomyPatientAware = corrected_asv_meta_stage.asv_meta_corrected.map { it }
            asvMetaForLungStatus = corrected_asv_meta_stage.asv_meta_corrected.map { it }
            asvMetaForMasterSummary = corrected_asv_meta_stage.asv_meta_corrected.map { it }
        }
    }
    if( bubbleplotterEnabled ) {
        BUBBLEPLOTTER(asvMetaForBubbleplotter)
    }
    if( umapClusteringEnabled ) {
        UMAP_CLUSTERING(asvMetaForUmap)
    }
    if( outlierEnabled ) {
        OUTLIER_CHECKER(
            asvClrForOutlier,
            metaMicroForOutlier
        )
    }
    if( collectorsEnabled ) {
        COLLECTORS_CURVE(
            asvFinalForCollectors,
            metaMicroForCollectors
        )
    }

    indicspecies_stage = null
    if( indicspeciesRunEarly ) {
        indicspecies_stage = INDICSPECIES(
            metaMicroForIndicspecies,
            asvFinalForIndicspecies
        )
        if( indicspeciesPlotEnabled ) {
            INDICSPECIES_PLOTS(
                metaMicroForIndicspeciesPlots,
                indicspecies_stage.all_tables.collect()
            )
        }
        if( indicspeciesAlignedEnabled ) {
            INDICSPECIES_ALIGNED_PLOTS(
                indicspecies_stage.all_tables.collect()
            )
        }
    }

    if( indicspeciesRunEarly ) {
        indicspeciesIsaForClustermaps = indicspecies_stage.all_tables.collect().map { tables ->
            def selected = clustermapsIsaAutoCandidates
                .collect { candidate -> tables.find { table -> table.name == candidate } }
                .find { it != null }
            if( selected == null ) {
                selected = tables
                    .findAll { table -> table.name.endsWith('_indicator_species_summary.tsv') }
                    .sort { table -> table.name }
                    .find()
            }
            if( selected == null ) {
                throw new IllegalStateException(
                    "INDICSPECIES produced no summary table suitable for CLUSTERMAPS. " +
                    "Expected one of: ${clustermapsIsaAutoCandidates.join(', ')}"
                )
            }
            selected
        }
    } else if( clustermapsExternalIsaTable ) {
        indicspeciesIsaForClustermaps = Channel.value(file(clustermapsExternalIsaTable))
    } else {
        indicspeciesIsaForClustermaps = Channel.value(file(emptyModulesPath))
    }
    indicspeciesReadyForPowerAnalysis = indicspeciesRunEarly ? indicspecies_stage.done.map { true } : Channel.value(false)
    indicspeciesTablesForOverlay = indicspeciesRunEarly ? indicspecies_stage.all_tables.collect() : Channel.value(file(emptyModulesPath))
    indicspeciesTablesForVocCorrelation = indicspeciesRunEarly ? indicspecies_stage.all_tables.collect() : Channel.value(file(emptyModulesPath))

    if( vocCorrelationEnabled ) {
        VOC_CORRELATION(
            asvMetaForVocCorrelation,
            asvFinalForVocCorrelation,
            indicspeciesTablesForVocCorrelation
        )
    }

    if( communityPredictorEnabled ) {
        COMMUNITY_PREDICTOR_COMPARISON(
            metaMicroForCommunityPredictor,
            asvFinalForCommunityPredictor,
            communityPredictorScriptHash
        )
    }

    if( diversityEnabled ) {
        def cohortInput = []
        if( diversityMatchedCohortSource == 'community_predictor_comparison' ) {
            if( !communityPredictorEnabled ) {
                exit 1, "Generated diversity matched cohort requires community_predictor_comparison.enabled=true"
            }
            cohortInput = COMMUNITY_PREDICTOR_COMPARISON.out.sample_cohort.ifEmpty {
                throw new IllegalStateException("Community predictor comparison produced no sample_analysis_cohort.tsv; cannot run matched-cohort diversity statistics")
            }
        } else if( diversityMatchedCohortPath ) {
            cohortInput = file(diversityMatchedCohortPath, checkIfExists: true)
        }
        DIVERSITY_ANALYSIS(metaMicroForDiversity, asvFinalForDiversity, cohortInput)
    }

    if( clustermapsEnabled ) {
        CLUSTERMAPS(
            asvMetaForClustermaps,
            metaMicroForClustermaps,
            indicspeciesIsaForClustermaps
        )
    }
    if( powerAnalysisEnabled ) {
        GROUP_POWER_ANALYSIS(
            asvMetaForPowerAnalysis,
            asvFinalForPowerAnalysis,
            indicspeciesReadyForPowerAnalysis
        )
    }
    if( taxonomyPatientAwareEnabled ) {
        TAXONOMY_GROUP_ASSOCIATION(
            asvMetaForTaxonomyPatientAware,
            asvFinalForTaxonomyPatientAware
        )
    }
    if( lungStatusAnalysisEnabled ) {
        PAIRED_GROUP_CONTRAST(
            asvMetaForLungStatus,
            asvFinalForLungStatus
        )
    }

    emit:
    meta_micro_network = metaMicroForNetwork
    asv_final_spieceasi = asvFinalForSpieceasi
    asv_final_network = asvFinalForNetwork
    asv_meta_master_summary = asvMetaForMasterSummary
    asv_final_master_summary = asvFinalForMasterSummary
    meta_micro_measurement_association = metaMicroForMeasurementAssociation
    asv_meta_measurement_association = asvMetaForMeasurementAssociation
    asv_final_measurement_association = asvFinalForMeasurementAssociation
    indicspecies_tables = indicspeciesTablesForOverlay
    indicspecies_group1_summary = indicspeciesRunEarly ? indicspecies_stage.group1_summary : Channel.value(file(emptyModulesPath))
}

workflow RUN_GRAPH_NETWORK {
    take:
    graph_all
    graph_thr
    node_features
    asv_final
    network_metadata
    asv_mag_done
    taxonomy_table
    isa_tables
    modules_sub
    modules_all
    microbial_state_done

    main:
    stage = GRAPH_NETWORK(
        graph_all,
        graph_thr,
        node_features,
        asv_final,
        network_metadata,
        asv_mag_done,
        taxonomy_table,
        isa_tables,
        modules_sub,
        modules_all,
        microbial_state_done
    )

    emit:
    done = stage.done
    module_best_stats_all = stage.module_best_stats_all
    layout_all = stage.layout_all
}

workflow RUN_ASV_MAG_NETWORK {
    take:
    graph
    node_features
    ecological_modules
    ecological_module_selection
    taxonomy_table
    asv_counts
    cruise_metadata
    sample_metadata
    isa_tables
    dep_asv_mag

    main:
    stage = ASV_MAG_NETWORK(
        graph,
        node_features,
        ecological_modules,
        ecological_module_selection,
        taxonomy_table,
        asv_counts,
        cruise_metadata,
        sample_metadata,
        isa_tables,
        dep_asv_mag
    )

    emit:
    done = stage.done
    analysis_mappings = stage.analysis_mappings
    functional_modules = stage.functional_modules
}

workflow RUN_MODULE_MAG_ANCHORS {
    take:
    modules_all
    node_features
    taxonomy_table
    asv_counts
    metadata_table
    asv_mag_mappings
    dep_asv_mag
    dep_graph_network

    main:
    stage = MODULE_MAG_ANCHORS(
        modules_all,
        node_features,
        taxonomy_table,
        asv_counts,
        metadata_table,
        asv_mag_mappings,
        dep_asv_mag,
        dep_graph_network
    )

    emit:
    asv_anchor_table = stage.asv_anchor_table
    module_summary = stage.module_summary
    sample_module_scores = stage.sample_module_scores
    sample_top_modules = stage.sample_top_modules
    sample_module_matrix = stage.sample_module_matrix
    done = stage.done
}

workflow RUN_GROUP_GUILD_FUNCTION {
    take:
    metadata_table
    asv_counts
    modules_all
    node_features
    isa_tables
    mappings
    functional_modules
    pca_scores
    pca_explained
    hybrid_assignments
    hybrid_centroids
    taxonomy_table

    main:
    stage = GROUP_GUILD_FUNCTION(
        metadata_table,
        asv_counts,
        modules_all,
        node_features,
        isa_tables,
        mappings,
        functional_modules,
        pca_scores,
        pca_explained,
        hybrid_assignments,
        hybrid_centroids,
        taxonomy_table,
        groupGuildFunctionScriptHash
    )

    emit:
    done = stage.done
}

def downloadReference(String downloadUrl, File destination) {
    destination.parentFile?.mkdirs()
    def tmpFile = File.createTempFile("sina_ref", ".download", destination.parentFile ?: new File('.'))
    tmpFile.withOutputStream { out ->
        new java.net.URL(downloadUrl).withInputStream { ins ->
            out << ins
        }
    }
    if( downloadUrl?.toLowerCase()?.endsWith('.gz') ) {
        destination.withOutputStream { out ->
            tmpFile.withInputStream { tmpIn ->
                new java.util.zip.GZIPInputStream(tmpIn).withCloseable { gz ->
                    out << gz
                }
            }
        }
        tmpFile.delete()
    } else {
        if( !tmpFile.renameTo(destination) ) {
            tmpFile.withInputStream { ins ->
                destination.withOutputStream { out ->
                    out << ins
                }
            }
            tmpFile.delete()
        }
    }
}

def ensureBlastReferenceExists(String basePath, String fastaPath, String label){
    if( fastaPath ) {
        def fasta = new File(fastaPath)
        if( !fasta.exists() ) {
            exit 1, "${label} FASTA not found: ${fastaPath}"
        }
        return
    }
    def baseFile = new File(basePath)
    if( baseFile.exists() || ['.nhr','.nin','.nsq','.ndb'].any { new File(basePath + it).exists() } ) {
        return
    }
    exit 1, "${label} BLAST database not found: ${basePath}. Provide a database prefix or a FASTA input."
}

def writeNormalizedManifest(List records, File destination, String sourceManifest){
    destination.parentFile?.mkdirs()
    def seen = new LinkedHashSet<String>()
    def lines = ['sample_id\tfastq_r1\tfastq_r2']
    records.each { rec ->
        def sampleId = rec.sample_id?.toString()?.trim()
        if( !sampleId || !seen.add(sampleId) ) {
            exit 1, "Duplicate or empty sample ID in ${sourceManifest ?: 'FASTQ discovery'}: ${sampleId}"
        }
        lines << [sampleId, new File(rec.r1.toString()).canonicalPath,
                  rec.paired && rec.r2 ? new File(rec.r2.toString()).canonicalPath : ''].join('\t')
    }
    destination.text = lines.join(System.lineSeparator()) + System.lineSeparator()
    return destination.canonicalPath
}

def safeFilename(String value){
    return value.replaceAll(/[^A-Za-z0-9._-]+/, '_').replaceAll(/^_+|_+$/, '') ?: 'group'
}

def generatedPaletteColor(int index, int total){
    // Golden-angle hue spacing remains stable when metadata row order is unchanged.
    float hue = ((index * 0.61803398875d) % 1.0d) as float
    int rgb = java.awt.Color.HSBtoRGB(hue, 0.62f, 0.78f)
    return String.format('#%06X', rgb & 0xFFFFFF)
}

def prepareMetadataAssets(String metadataPath, String sampleCol, String groupCol, String colorCol,
                          String palettePath, File outputBase){
    File source = new File(metadataPath)
    def rows = source.readLines('UTF-8').findAll { it != null && !it.trim().isEmpty() }
    if( !rows ) {
        exit 1, "Metadata table is empty: ${metadataPath}"
    }
    def header = rows[0].split(/\t/, -1).collect { it.trim() }
    int sampleIdx = header.indexOf(sampleCol)
    int groupIdx = header.indexOf(groupCol)
    if( sampleIdx < 0 || groupIdx < 0 ) {
        exit 1, "Metadata must contain configured columns '${sampleCol}' and '${groupCol}': ${metadataPath}"
    }
    int colorIdx = header.indexOf(colorCol)
    def groups = []
    rows.drop(1).each { line ->
        def fields = line.split(/\t/, -1)
        if( fields.length > groupIdx ) {
            def value = fields[groupIdx].trim()
            if( value && !groups.contains(value) ) groups << value
        }
    }
    LinkedHashMap<String,String> palette = [:]
    if( palettePath ) {
        File paletteFile = new File(palettePath)
        if( !paletteFile.exists() ) exit 1, "Metadata palette file not found: ${palettePath}"
        def paletteRows = paletteFile.readLines('UTF-8').findAll { it?.trim() }
        paletteRows.eachWithIndex { line, idx ->
            def fields = line.split(/\t|,/, -1).collect { it.trim() }
            if( fields.size() >= 2 && !(idx == 0 && fields[0].equalsIgnoreCase('value')) ) {
                palette[fields[0]] = fields[1]
            }
        }
    }
    groups.eachWithIndex { group, idx ->
        if( !palette[group] && colorIdx >= 0 ) {
            def matching = rows.drop(1).find { line ->
                def fields = line.split(/\t/, -1)
                fields.length > Math.max(groupIdx, colorIdx) && fields[groupIdx].trim() == group && fields[colorIdx].trim()
            }
            if( matching ) palette[group] = matching.split(/\t/, -1)[colorIdx].trim()
        }
        if( !palette[group] ) palette[group] = generatedPaletteColor(idx, groups.size())
    }
    // Palette files may be shared across studies; publish only observed groups.
    LinkedHashMap<String,String> observedPalette = [:]
    groups.each { group -> observedPalette[group] = palette[group] }
    palette = observedPalette
    outputBase.parentFile?.mkdirs()
    File paletteOut = new File(outputBase.parentFile, "${outputBase.name}_palette.tsv")
    paletteOut.text = 'value\tcolor' + System.lineSeparator() + palette.collect { key, value -> "${key}\t${value}" }.join(System.lineSeparator()) + System.lineSeparator()
    File metadataOut = new File(outputBase.parentFile, "${outputBase.name}_metadata.tsv")
    def outputHeader = colorIdx >= 0 ? header : header + [colorCol]
    def outputRows = [outputHeader.join('\t')]
    rows.drop(1).each { line ->
        def fields = line.split(/\t/, -1).toList()
        def missingFields = header.size() - fields.size()
        if( missingFields > 0 ) {
            fields.addAll((1..missingFields).collect { '' })
        }
        def group = fields[groupIdx].trim()
        if( colorIdx >= 0 ) fields[colorIdx] = palette[group] ?: fields[colorIdx]
        else fields << (palette[group] ?: '')
        outputRows << fields.join('\t')
    }
    metadataOut.text = outputRows.join(System.lineSeparator()) + System.lineSeparator()
    return [metadata: metadataOut.canonicalPath, palette: paletteOut.canonicalPath]
}

def shellQuote(String value){
    if( value == null ){
        return "''"
    }
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def joinShellArgs(List paths){
    if( !paths ) {
        return ''
    }
    return paths.collect { shellQuote(it.toString()) }.join(' ')
}

/**
 * Helpers
 */
def normalizeList(value, fallback){
    if( !value ) return fallback
    if( value instanceof List ) return value.collect { it.toString() }
    return value.toString().split(/\|/).collect { it.trim() }.findAll { it }
}

def normalizePresetList(value, fallback=[], presets=[:]){
    if( !value ) return fallback
    if( value instanceof List ) {
        return value.collect { it.toString().trim() }.findAll { it }
    }

    def text = value.toString().trim()
    if( !text ) return fallback

    def presetValue = presets[text]
    if( presetValue == null && text.endsWith('_order') ) {
        presetValue = presets[text.replaceFirst(/_order$/, '')]
    }
    if( presetValue != null && presetValue != value ) {
        return normalizePresetList(presetValue, fallback, presets)
    }

    return text.split(/[,|]/).collect { it.trim() }.findAll { it }
}

def compilePatterns(value, fallback){
    def list = value ?: fallback
    return list.collect { java.util.regex.Pattern.compile(it.toString()) }
}

def matchesExtension(String name, List<java.util.regex.Pattern> patterns){
    patterns.any { it.matcher(name).find() }
}

def isR1Like(String base, List<String> tokens){
    tokens.any { tok ->
        def rx = java.util.regex.Pattern.compile("(^|[_\\.\\-])${java.util.regex.Pattern.quote(tok)}([_\\.\\-]|\$)")
        rx.matcher(base).find()
    }
}

def sampleFromName(String baseName, String stripRegex, List<java.util.regex.Pattern> extPatterns){
    def base = baseName
    extPatterns.each { base = base.replaceAll(it, '') }
    base = base.replaceAll(stripRegex, '')
    base = base.replaceAll(/[_\-.]+$/, '')
    return base
}

def findR2File(File r1File, List<String> r1Tokens, List<String> r2Tokens){
    def original = r1File.name
    File matched = null
    (0..<r1Tokens.size()).each { i ->
        if( matched != null ) {
            return
        }
        def r1 = r1Tokens[i]
        def r2 = r2Tokens[i]
        def replacements = [
            ["_${java.util.regex.Pattern.quote(r1)}_", "_${r2}_"],
            ["\\.${java.util.regex.Pattern.quote(r1)}\\.", ".${r2}."],
            ["-${java.util.regex.Pattern.quote(r1)}-", "-${r2}-"],
            ["-${java.util.regex.Pattern.quote(r1)}\\.", "-${r2}."],
            ["_${java.util.regex.Pattern.quote(r1)}\\.", "_${r2}."],
            ["_${java.util.regex.Pattern.quote(r1)}\$", "_${r2}"],
            ["${java.util.regex.Pattern.quote(r1)}_001", "${r2}_001"]
        ]
        replacements.each { rep ->
            if( matched != null ) {
                return
            }
            def candidateName = original.replaceFirst(rep[0], rep[1])
            if( candidateName != original ){
                def candidate = new File(r1File.parentFile, candidateName)
                if( candidate.exists() ) {
                    matched = candidate
                }
            }
        }
    }
    return matched
}

def collectSampleRecords(String inputDirPath, List<String> r1Tokens, List<String> r2Tokens,
                         List<java.util.regex.Pattern> extPatterns, String stripRegex, boolean allowSingleEnd){
    File dir = new File(inputDirPath)
    if( !dir.exists() ){
        throw new IllegalArgumentException("Input directory does not exist: ${inputDirPath}")
    }
    def files = dir.listFiles()?.findAll { it.isFile() && matchesExtension(it.name, extPatterns) }?.sort { a, b -> a.name <=> b.name } ?: []
    def records = []
    files.each { file ->
        if( isR1Like(file.name, r1Tokens) ){
            def sampleId = sampleFromName(file.name, stripRegex, extPatterns)
            def r2File = findR2File(file, r1Tokens, r2Tokens)
            if( r2File ){
                records << [ sample_id: sampleId, paired: true, r1: file.canonicalPath, r2: r2File.canonicalPath ]
            } else if( allowSingleEnd ) {
                records << [ sample_id: sampleId, paired: false, r1: file.canonicalPath ]
            } else {
                log.warn "Skipping ${file.name}: no matching R2 detected."
            }
        }
    }
    return records
}

def loadManifestSamples(String manifestPath){
    File manifest = new File(manifestPath)
    if( !manifest.exists() ) {
        exit 1, "Manifest file not found: ${manifestPath}"
    }
    def records = []
    manifest.eachLine { line ->
        def trimmed = line.trim()
        if( !trimmed || trimmed.startsWith('#') ) {
            return
        }
        def parts = trimmed.split(/\t/)
        if( parts.length >= 2 && parts[0].trim().equalsIgnoreCase('sample_id') &&
            parts[1].trim().toLowerCase() in ['fastq_r1', 'r1', 'read1'] ) {
            return
        }
        if( parts.length < 2 ) {
            exit 1, "Manifest line must contain sample_id and R1 path separated by tab: ${line}"
        }
        def sampleId = parts[0].trim()
        def r1Path = resolveOptionalPath(parts[1].trim(), manifest.parentFile)
        def r2Path = parts.length > 2 && parts[2].trim() ? resolveOptionalPath(parts[2].trim(), manifest.parentFile) : null
        if( !sampleId ) {
            exit 1, "Sample ID missing in manifest line: ${line}"
        }
        if( !r1Path ) {
            exit 1, "R1 path missing in manifest line: ${line}"
        }
        def r1File = new File(r1Path)
        if( !r1File.exists() ) {
            exit 1, "R1 file not found for sample ${sampleId}: ${r1Path}"
        }
        File r2File = null
        boolean paired = false
        if( r2Path ) {
            r2File = new File(r2Path)
            if( !r2File.exists() ) {
                exit 1, "R2 file not found for sample ${sampleId}: ${r2Path}"
            }
            paired = true
        }
        records << [
            sample_id: sampleId,
            paired: paired,
            r1: r1File.canonicalPath,
            r2: paired ? r2File.canonicalPath : null
        ]
    }
    return records
}

def resolveOptionalPath(String pathValue, File baseDir){
    if( !pathValue ) return null
    def candidate = new File(pathValue)
    if( candidate.isAbsolute() ) {
        return candidate.canonicalPath
    }
    def baseCandidate = new File(baseDir ?: new File('.'), pathValue)
    if( baseCandidate.exists() ) {
        return baseCandidate.canonicalPath
    }
    def projectCandidate = new File(new File(projectDir.toString()), pathValue)
    if( projectCandidate.exists() ) {
        return projectCandidate.canonicalPath
    }
    return baseCandidate.canonicalPath
}

def resolveOutputRelative(String pathValue, String baseOutputDir){
    if( !pathValue ) return baseOutputDir
    def candidate = new File(pathValue)
    if( candidate.isAbsolute() ) {
        return candidate.canonicalPath
    }
    return new File(baseOutputDir ?: '.', pathValue).canonicalPath
}

def extractNamedStringMap(Map cfg, List<String> orderedKeys, String nestedKey, String legacySuffix){
    LinkedHashMap<String,String> out = [:]
    if( cfg?.get(nestedKey) instanceof Map ) {
        cfg[nestedKey].each { k, v ->
            def key = k?.toString()?.trim()
            def val = v?.toString()?.trim()
            if( key && val ) {
                out[key] = val
            }
        }
    }
    orderedKeys.eachWithIndex { key, idx ->
        def legacyKey = "group${idx + 1}_${legacySuffix}"
        if( cfg?.containsKey(legacyKey) ) {
            def val = cfg[legacyKey]?.toString()?.trim()
            if( val ) {
                out[key] = val
            }
        }
    }
    return out
}

def extractNamedListMap(Map cfg, List<String> orderedKeys, String nestedKey, String legacySuffix){
    LinkedHashMap<String,List<String>> out = [:]
    if( cfg?.get(nestedKey) instanceof Map ) {
        cfg[nestedKey].each { k, v ->
            def key = k?.toString()?.trim()
            def vals = []
            if( v instanceof List ) {
                vals = v.collect { it?.toString()?.trim() }.findAll { it }
            } else if( v ) {
                vals = v.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
            }
            if( key && vals ) {
                out[key] = vals
            }
        }
    }
    orderedKeys.eachWithIndex { key, idx ->
        def legacyKey = "group${idx + 1}_${legacySuffix}"
        if( cfg?.containsKey(legacyKey) ) {
            def raw = cfg[legacyKey]
            def vals = []
            if( raw instanceof List ) {
                vals = raw.collect { it?.toString()?.trim() }.findAll { it }
            } else if( raw ) {
                vals = raw.toString().split(/[,|]/).collect { it.trim() }.findAll { it }
            }
            if( vals ) {
                out[key] = vals
            }
        }
    }
    return out
}
