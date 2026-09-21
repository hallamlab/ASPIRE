#!/usr/bin/env bash
set -euo pipefail

export NXF_VER="${NXF_VER:-25.10.0}"
export NXF_SYNTAX_PARSER="${NXF_SYNTAX_PARSER:-v1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "$0")"
CONTROL_ENV_DIR="${CONTROL_ENV_DIR:-${SCRIPT_DIR}/.controller_env}"
ORIGINAL_ARGS=("$@")

PROCESS_ORDER=(
  FASTP_QC
  MERGE_READS
  FILTER_READS
  RELABEL_FILTERED
  CONCAT_FASTAS
  DEREPLICATE
  DENOISE
  CHIMERA_CHECK
  CREATE_COUNT_MATRIX
  FILTER_TABLE
  SINA_TRIM
  TAXONOMY
  PREPARE_BLAST_DATABASES
  MITOMASTER
  MITO_DECONTAM
  FILTER_COUNTS
  GENERAL_STATS
  PREPROCESS_DATASET
  VALIDATE_PREPROCESS_DATASET
  PLOT_METADATA
  ASV_TIME_DEPTH_CURTAIN
  GROUPING_DIAGNOSTICS
  GROUP_LABEL_AUGMENTATION
  PLOT_UPSET
  ASV_BATCH_CORRECTION
  ASV_META_FROM_CORRECTED
  BUBBLEPLOTTER
  UMAP_CLUSTERING
  OUTLIER_CHECKER
  COLLECTORS_CURVE
  DIVERSITY_ANALYSIS
  INDICSPECIES
  INDICSPECIES_PLOTS
  INDICSPECIES_ALIGNED_PLOTS
  VOC_CORRELATION
  COMMUNITY_PREDICTOR_COMPARISON
  CLUSTERMAPS
  GROUP_POWER_ANALYSIS
  TAXONOMY_GROUP_ASSOCIATION
  PAIRED_GROUP_CONTRAST
  ASV_MAG_LINK
  SPIECEASI
  MEASUREMENT_ASSOCIATION
  TITAN_INSTALL
  TITAN_PREPARE
  TITAN_ANALYSIS
  TITAN_COLLECT
  TITAN_PLOTS
  MICROBIAL_COMPARTMENT_PREPARE
  MICROBIAL_COMPARTMENT_INFERENCE
  MICROBIAL_COMPARTMENT_POSTHOC
  MICROBIAL_COMPARTMENT_PLOTS
  COMMUNITY_TURNOVER_PREPARE
  COMMUNITY_TURNOVER_ANALYSIS
  COMMUNITY_TURNOVER_LCBD
  COMMUNITY_TURNOVER_PLOTS
  NETWORK_MODULES
  MICROBIAL_STATE_INTERPRETATION
  GRAPH_NETWORK
  MODULE_MEASUREMENT_ASSOCIATION
  ASV_MAG_NETWORK
  GENOME_COOCCURRENCE
  ASV_MAG_CURTAINS
  MODULE_MAG_ANCHORS
  GROUP_GUILD_FUNCTION
  ECOLOGICAL_CONTEXT_ATLAS
  SANKEY
  MASTER_SUMMARY
)

declare -A PROCESS_ALIASES=(
  [POWER_ANALYSIS_PIPELINE]=GROUP_POWER_ANALYSIS
  [TAXONOMY_PATIENT_AWARE]=TAXONOMY_GROUP_ASSOCIATION
  [LUNG_STATUS_ANALYSIS]=PAIRED_GROUP_CONTRAST
)

usage() {
  cat <<'EOF'
Usage:
  run_asv_pipeline.sh [CONFIG_FILE] [--phase PHASE] [--rerun-from PROCESS_NAME] [--resume-run RUN_NAME] [--no-resume] [--list-stages] [-- NEXTFLOW_ARGS...]
  run_asv_pipeline.sh [CONFIG_FILE] cache list [--limit N]
  run_asv_pipeline.sh [CONFIG_FILE] cache usage
  run_asv_pipeline.sh [CONFIG_FILE] cache prune [--keep N] [--force]
  run_asv_pipeline.sh [CONFIG_FILE] cache clear (--run RUN_NAME | --all) [--include-conda] [--force]

Options:
  --phase PHASE             Run all, preprocess, or analysis (default: all).
                            Analysis consumes the canonical preprocessing dataset
                            and never instantiates raw-read processing.
  --rerun-from PROCESS_NAME  Resume the normal workflow while forcing only PROCESS_NAME
                             to execute. True dependency descendants rerun automatically;
                             independent branches remain cacheable.
  --resume-run RUN_NAME      Resume from a specific Nextflow run name/id (from `nextflow log -q`).
  --no-resume                Disable resume for this run (cold execution).
  --prune-cache-after-run N  After success, retain N run records and clean older,
                             unreferenced Nextflow work (recommended: 10).
  --list-stages              Print known process names for --rerun-from and exit.
  --help, -h                 Show this help.

Cache management:
  cache list                 List recent Nextflow run records and protected phase pointers.
  cache usage                Report disk use for work, Conda, Nextflow metadata, and staging.
  cache prune --keep N       Retain N run records, including protected phase pointers.
                             Dry-run unless --force is supplied.
  cache clear --run NAME     Remove one run record/cache entry; dry-run unless --force.
  cache clear --all          Select every run record/cache entry; dry-run unless --force.
  --include-conda            With `cache clear --all --force`, also remove shared Conda
                             environments/packages so they are rebuilt on the next run.

Examples:
  run_asv_pipeline.sh asv_pipeline_nextflow.yml --phase preprocess
  run_asv_pipeline.sh asv_pipeline_nextflow.yml --phase analysis
  run_asv_pipeline.sh asv_pipeline_nextflow.yml cache prune --keep 10
  run_asv_pipeline.sh asv_pipeline_nextflow.yml cache prune --keep 10 --force
  run_asv_pipeline.sh asv_pipeline_nextflow.yml --rerun-from FILTER_COUNTS
  run_asv_pipeline.sh --resume-run lethal_poisson
  run_asv_pipeline.sh --rerun-from PLOT_METADATA -- -with-report report.html
EOF
}

if [[ -z "${IN_CONTROLLER_ENV:-}" ]]; then
  if ! command -v mamba >/dev/null 2>&1; then
    echo "mamba is required to bootstrap the controller environment." >&2
    exit 1
  fi
  exec 8>"${TMPDIR:-/tmp}/aspire-controller-${UID}.lock"
  if ! flock -w 7200 8; then
    echo "Timed out waiting for another ASPIRE controller environment update." >&2
    exit 1
  fi
  if [[ ! -d "$CONTROL_ENV_DIR" ]]; then
    echo "[controller] Creating mamba env at $CONTROL_ENV_DIR"
    mamba env create --yes --prefix "$CONTROL_ENV_DIR" --file "${SCRIPT_DIR}/processes/controller/env.yml"
  elif [[ ! -x "$CONTROL_ENV_DIR/bin/mamba" ]]; then
    echo "[controller] Updating controller environment to provide an isolated mamba executable."
    mamba env update --prefix "$CONTROL_ENV_DIR" --file "${SCRIPT_DIR}/processes/controller/env.yml"
  fi
  flock -u 8
  exec env IN_CONTROLLER_ENV=1 CONTROL_ENV_DIR="$CONTROL_ENV_DIR" conda run --no-capture-output -p "$CONTROL_ENV_DIR" "$SCRIPT_PATH" "$@"
fi

CONFIG_FILE="asv_pipeline_nextflow.yml"
CONFIG_SET=0
LIST_STAGES=0
RERUN_FROM=""
NEXTFLOW_ARGS=()
RESUME_ENABLED=1
RESUME_RUN=""
ASPIRE_PHASE="all"
CACHE_MODE=0
CACHE_ACTION=""
CACHE_KEEP=10
CACHE_LIMIT=20
CACHE_FORCE=0
CACHE_CLEAR_RUN=""
CACHE_CLEAR_ALL=0
CACHE_INCLUDE_CONDA=0
PRUNE_CACHE_AFTER_RUN=0

if [[ "${1:-}" == "cache" || ( $# -ge 2 && "${2:-}" == "cache" && "${1:-}" != -* ) ]]; then
  CACHE_MODE=1
  if [[ "${1:-}" == "cache" ]]; then
    shift
  else
    CONFIG_FILE="$1"
    CONFIG_SET=1
    shift 2
  fi
  CACHE_ACTION="${1:-list}"
  [[ $# -gt 0 ]] && shift
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --keep)
        [[ $# -ge 2 ]] || { echo "--keep requires a positive integer" >&2; exit 1; }
        CACHE_KEEP="$2"
        shift 2
        ;;
      --keep=*) CACHE_KEEP="${1#*=}"; shift ;;
      --limit)
        [[ $# -ge 2 ]] || { echo "--limit requires a positive integer" >&2; exit 1; }
        CACHE_LIMIT="$2"
        shift 2
        ;;
      --limit=*) CACHE_LIMIT="${1#*=}"; shift ;;
      --run)
        [[ $# -ge 2 ]] || { echo "--run requires a Nextflow run name" >&2; exit 1; }
        CACHE_CLEAR_RUN="$2"
        shift 2
        ;;
      --run=*) CACHE_CLEAR_RUN="${1#*=}"; shift ;;
      --all) CACHE_CLEAR_ALL=1; shift ;;
      --include-conda) CACHE_INCLUDE_CONDA=1; shift ;;
      --force) CACHE_FORCE=1; shift ;;
      -h|--help) usage; exit 0 ;;
      *) echo "Unknown cache option: $1" >&2; exit 1 ;;
    esac
  done
fi

while [[ "$CACHE_MODE" -eq 0 && $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --list-stages)
      LIST_STAGES=1
      shift
      ;;
    --phase)
      if [[ $# -lt 2 ]]; then
        echo "--phase requires one of: all, preprocess, analysis" >&2
        exit 1
      fi
      ASPIRE_PHASE="${2,,}"
      shift 2
      ;;
    --phase=*)
      ASPIRE_PHASE="${1#*=}"
      ASPIRE_PHASE="${ASPIRE_PHASE,,}"
      shift
      ;;
    --rerun-from)
      if [[ $# -lt 2 ]]; then
        echo "--rerun-from requires a process name" >&2
        exit 1
      fi
      RERUN_FROM="$2"
      shift 2
      ;;
    --rerun-from=*)
      RERUN_FROM="${1#*=}"
      shift
      ;;
    --resume-run)
      if [[ $# -lt 2 ]]; then
        echo "--resume-run requires a run name or id" >&2
        exit 1
      fi
      RESUME_RUN="$2"
      shift 2
      ;;
    --resume-run=*)
      RESUME_RUN="${1#*=}"
      shift
      ;;
    --no-resume)
      RESUME_ENABLED=0
      shift
      ;;
    --prune-cache-after-run)
      if [[ $# -lt 2 || ! "$2" =~ ^[1-9][0-9]*$ ]]; then
        echo "--prune-cache-after-run requires a positive integer" >&2
        exit 1
      fi
      PRUNE_CACHE_AFTER_RUN="$2"
      shift 2
      ;;
    --prune-cache-after-run=*)
      PRUNE_CACHE_AFTER_RUN="${1#*=}"
      [[ "$PRUNE_CACHE_AFTER_RUN" =~ ^[1-9][0-9]*$ ]] || { echo "--prune-cache-after-run requires a positive integer" >&2; exit 1; }
      shift
      ;;
    --)
      shift
      NEXTFLOW_ARGS+=("$@")
      break
      ;;
    -*)
      NEXTFLOW_ARGS+=("$1")
      shift
      ;;
    *)
      if [[ $CONFIG_SET -eq 0 ]]; then
        CONFIG_FILE="$1"
        CONFIG_SET=1
      else
        NEXTFLOW_ARGS+=("$1")
      fi
      shift
  esac
done

if [[ "$CACHE_MODE" -eq 1 ]]; then
  case "$CACHE_ACTION" in
    list|status|usage) ;;
    prune)
      [[ "$CACHE_KEEP" =~ ^[1-9][0-9]*$ ]] || { echo "--keep must be a positive integer" >&2; exit 1; }
      ;;
    clear)
      if [[ "$CACHE_CLEAR_ALL" -eq 1 && -n "$CACHE_CLEAR_RUN" ]]; then
        echo "cache clear accepts either --run or --all, not both" >&2
        exit 1
      fi
      if [[ "$CACHE_CLEAR_ALL" -eq 0 && -z "$CACHE_CLEAR_RUN" ]]; then
        echo "cache clear requires --run RUN_NAME or --all" >&2
        exit 1
      fi
      if [[ "$CACHE_INCLUDE_CONDA" -eq 1 && "$CACHE_CLEAR_ALL" -ne 1 ]]; then
        echo "--include-conda is valid only with cache clear --all" >&2
        exit 1
      fi
      ;;
    *) echo "Unknown cache command: ${CACHE_ACTION}" >&2; exit 1 ;;
  esac
  [[ "$CACHE_LIMIT" =~ ^[1-9][0-9]*$ ]] || { echo "--limit must be a positive integer" >&2; exit 1; }
fi

case "$ASPIRE_PHASE" in
  all|preprocess|analysis) ;;
  *)
    echo "--phase requires one of: all, preprocess, analysis (received: ${ASPIRE_PHASE})" >&2
    exit 1
    ;;
esac

if [[ "$LIST_STAGES" -eq 1 ]]; then
  printf '%s\n' "${PROCESS_ORDER[@]}"
  if [[ ${#PROCESS_ALIASES[@]} -gt 0 ]]; then
    echo "# aliases"
    for alias_name in "${!PROCESS_ALIASES[@]}"; do
      printf '%s -> %s\n' "$alias_name" "${PROCESS_ALIASES[$alias_name]}"
    done | sort
  fi
  exit 0
fi

if [[ "$RESUME_ENABLED" -eq 0 && -n "$RESUME_RUN" ]]; then
  echo "--no-resume cannot be combined with --resume-run" >&2
  exit 1
fi
if [[ "$RESUME_ENABLED" -eq 0 && -n "$RERUN_FROM" ]]; then
  echo "--no-resume cannot be combined with --rerun-from; stage-aware reruns require the normal Nextflow resume cache" >&2
  exit 1
fi
sanitize_nextflow_args() {
  local sanitized=()
  local skip_next=0
  local i arg next_arg
  for ((i=0; i<${#NEXTFLOW_ARGS[@]}; i++)); do
    if [[ $skip_next -eq 1 ]]; then
      skip_next=0
      continue
    fi
    arg="${NEXTFLOW_ARGS[$i]}"
    case "$arg" in
      -resume)
        if (( i + 1 < ${#NEXTFLOW_ARGS[@]} )); then
          next_arg="${NEXTFLOW_ARGS[$((i + 1))]}"
          if [[ "$next_arg" != -* ]]; then
            skip_next=1
          fi
        fi
        echo "[controller] Ignoring passthrough '-resume' argument; use --resume-run or --no-resume." >&2
        ;;
      -resume=*)
        echo "[controller] Ignoring passthrough '-resume=...' argument; use --resume-run or --no-resume." >&2
        ;;
      *)
        sanitized+=("$arg")
        ;;
    esac
  done
  NEXTFLOW_ARGS=("${sanitized[@]}")
}

successful_task_count() {
  local run_name="$1"
  local rows rc process_name workdir status count=0
  set +e
  rows="$(nextflow log "$run_name" -f 'process,workdir,status' 2>&1)"
  rc=$?
  set -e
  if [[ $rc -ne 0 ]]; then
    echo 0
    return 0
  fi
  while IFS=$'\t' read -r process_name workdir status; do
    [[ -n "$process_name" ]] || continue
    [[ "$process_name" == "process" && "$workdir" == "workdir" ]] && continue
    case "${status^^}" in
      COMPLETED|CACHED) count=$((count + 1)) ;;
    esac
  done <<< "$rows"
  echo "$count"
}

select_baseline_run() {
  local selected=""
  local phase_pointer="${RUNTIME_DIR}/last_successful_run.${ASPIRE_PHASE}"
  local runs=()
  local i run_name task_count best_count=-1

  if [[ -n "$RESUME_RUN" ]]; then
    if ! nextflow log "$RESUME_RUN" >/dev/null 2>&1; then
      echo "Specified --resume-run not found in Nextflow history: ${RESUME_RUN}" >&2
      exit 1
    fi
    selected="$RESUME_RUN"
    echo "$selected"
    return 0
  fi

  if [[ -f "$phase_pointer" ]]; then
    selected="$(head -n 1 "$phase_pointer" | tr -d '[:space:]')"
    if [[ -n "$selected" ]] && nextflow log "$selected" >/dev/null 2>&1; then
      echo "$selected"
      return 0
    fi
    echo "[controller] Ignoring stale phase resume pointer: ${phase_pointer}" >&2
  fi

  if [[ "$ASPIRE_PHASE" != "all" ]]; then
    echo ""
    return 0
  fi

  mapfile -t runs < <(cache_managed_runs)
  if [[ ${#runs[@]} -eq 0 ]]; then
    echo ""
    return 0
  fi
  # Prefer the most complete cache lineage, not merely the newest interrupted
  # attempt. Iterating newest-to-oldest without replacing ties selects the
  # newest run among those with the same successful-task count.
  for ((i=${#runs[@]}-1; i>=0; i--)); do
    run_name="${runs[$i]}"
    task_count="$(successful_task_count "$run_name")"
    if (( task_count > best_count )); then
      selected="$run_name"
      best_count=$task_count
    fi
  done

  if [[ -n "$selected" && $best_count -gt 0 ]]; then
    echo "$selected"
  else
    echo ""
  fi
}

sanitize_nextflow_args

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  exit 1
fi
CONFIG_FILE="$(realpath "$CONFIG_FILE")"
CONFIG_DIR="$(dirname "$CONFIG_FILE")"

GROUPING_REQUIRE_CRUISE_LEVEL=$(yq -r '.grouping_diagnostics.require_cruise_level_groups // false' "$CONFIG_FILE")
GROUPING_CRUISE_LEVEL_COUNT=$(yq -r '[.grouping_diagnostics.cruise_level_group_cols[]? | select(. != null and (tostring | length) > 0)] | length' "$CONFIG_FILE")
if [[ "$GROUPING_REQUIRE_CRUISE_LEVEL" == "true" && "$GROUPING_CRUISE_LEVEL_COUNT" -eq 0 ]]; then
  echo "grouping_diagnostics.require_cruise_level_groups=true requires at least one cruise_level_group_cols entry" >&2
  exit 1
fi

resolve_config_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    realpath -m "$value"
  else
    realpath -m "${CONFIG_DIR}/${value}"
  fi
}

OUTPUT_DIR=$(yq -r '.paths.output_dir // empty' "$CONFIG_FILE")
RUNTIME_DIR=$(yq -r '.paths.runtime_dir // empty' "$CONFIG_FILE")
KEEP_RUNTIME_DIR=$(yq -r '.paths.keep_runtime_dir // true' "$CONFIG_FILE")
WORK_DIR=$(yq -r '.paths.work_dir // empty' "$CONFIG_FILE")
CONDA_CACHE_DIR=$(yq -r '.paths.conda_cache_dir // empty' "$CONFIG_FILE")

if [[ -z "$OUTPUT_DIR" ]]; then
  echo "paths.output_dir must be set in $CONFIG_FILE" >&2
  exit 1
fi

OUTPUT_DIR="$(resolve_config_path "$OUTPUT_DIR")"
if [[ -z "$RUNTIME_DIR" || "$RUNTIME_DIR" == "null" ]]; then
  RUNTIME_DIR="${OUTPUT_DIR}/.aspire"
else
  RUNTIME_DIR="$(resolve_config_path "$RUNTIME_DIR")"
fi
if [[ -z "$WORK_DIR" || "$WORK_DIR" == "null" ]]; then
  WORK_DIR="${RUNTIME_DIR}/nf_work"
else
  WORK_DIR="$(resolve_config_path "$WORK_DIR")"
fi
if [[ -z "$CONDA_CACHE_DIR" || "$CONDA_CACHE_DIR" == "null" ]]; then
  CONDA_CACHE_DIR="${RUNTIME_DIR}/conda_cache"
else
  CONDA_CACHE_DIR="$(resolve_config_path "$CONDA_CACHE_DIR")"
fi
PUBLICATION_STAGING_DIR="${RUNTIME_DIR}/publication_staging"

case "${KEEP_RUNTIME_DIR,,}" in
  true|false) ;;
  *)
    echo "paths.keep_runtime_dir must be true or false: ${KEEP_RUNTIME_DIR}" >&2
    exit 1
    ;;
esac

mkdir -p "$WORK_DIR"
mkdir -p "$CONDA_CACHE_DIR"
mkdir -p "${CONDA_CACHE_DIR}/pkgs"
mkdir -p \
  "${OUTPUT_DIR}/modules" \
  "${OUTPUT_DIR}/intermediates" \
  "${OUTPUT_DIR}/references" \
  "${OUTPUT_DIR}/summary" \
  "${OUTPUT_DIR}/logs"

# Mutating cache commands and workflow executions require exclusive ownership.
# A read-only cache listing remains available while a workflow is active.
if [[ "$CACHE_MODE" -eq 0 || "$CACHE_ACTION" != "list" && "$CACHE_ACTION" != "status" && "$CACHE_ACTION" != "usage" ]]; then
  exec 7>"${CONDA_CACHE_DIR}/.aspire-run.lock"
  if ! flock -n 7; then
    echo "Another ASPIRE run is actively using Conda cache: ${CONDA_CACHE_DIR}" >&2
    echo "Wait for that run to finish before executing a workflow or mutating its cache." >&2
    exit 1
  fi
fi

cache_managed_runs() {
  # Nextflow history is repository-wide and can include tests or executions
  # using other output trees.  Cache management is intentionally limited to
  # runs that used the work directory resolved from the supplied config.
  nextflow log 2>/dev/null | awk -F '\t' -v expected_work="-w ${WORK_DIR}" '
    function trim(value) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
      return value
    }
    index($7, expected_work) > 0 {
      run_name = trim($3)
      if (run_name != "") print run_name
    }
  '
}

cache_run_is_managed() {
  local requested="$1" run_name
  while IFS= read -r run_name; do
    [[ "$run_name" == "$requested" ]] && return 0
  done < <(cache_managed_runs)
  return 1
}

cache_pointer_runs() {
  local pointer run_name
  for pointer in "${RUNTIME_DIR}"/last_successful_run.*; do
    [[ -f "$pointer" ]] || continue
    run_name="$(head -n 1 "$pointer" | tr -d '[:space:]')"
    [[ -n "$run_name" ]] && printf '%s\n' "$run_name"
  done | sort -u
}

legacy_successful_pipeline_run() {
  # Phase-aware pointers were introduced after existing ASPIRE histories had
  # already accumulated successful runs.  During that one-time migration,
  # identify only the newest successful, non-preview ASPIRE run that used this
  # configuration's work directory.  This excludes controller/unit-test runs
  # and prevents cache pruning from adopting an unrelated Nextflow lineage.
  nextflow log 2>/dev/null | awk -F '\t' -v expected_work="-w ${WORK_DIR}" '
    function trim(value) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
      return value
    }
    {
      run_name = trim($3)
      run_status = trim($4)
      command = $7
      if (run_status == "OK" &&
          index(command, expected_work) > 0 &&
          command !~ /(^|[[:space:]])-preview([[:space:]]|$)/ &&
          command ~ /asv_pipeline/) {
        selected = run_name
      }
    }
    END {
      if (selected != "") print selected
    }
  '
}

adopt_legacy_cache_baseline() {
  local selected pointer_tmp
  selected="$(legacy_successful_pipeline_run)"
  [[ -n "$selected" ]] || return 1
  nextflow log "$selected" >/dev/null 2>&1 || return 1

  pointer_tmp="${RUNTIME_DIR}/last_successful_run.all.tmp.$$"
  printf '%s\n' "$selected" > "$pointer_tmp"
  mv "$pointer_tmp" "${RUNTIME_DIR}/last_successful_run.all"
  echo "[cache] Adopted legacy successful run '${selected}' as the protected all-phase baseline." >&2
  return 0
}

cache_list_runs() {
  local runs=() protected=() start=0 i run_name marker
  mapfile -t runs < <(cache_managed_runs)
  mapfile -t protected < <(cache_pointer_runs)
  if (( ${#runs[@]} > CACHE_LIMIT )); then
    start=$((${#runs[@]} - CACHE_LIMIT))
  fi
  printf 'run_name\tprotected_phase_resume\n'
  for ((i=start; i<${#runs[@]}; i++)); do
    run_name="${runs[$i]}"
    marker="no"
    if printf '%s\n' "${protected[@]:-}" | grep -Fxq "$run_name"; then marker="yes"; fi
    printf '%s\t%s\n' "$run_name" "$marker"
  done
  printf '\nManaged run records for configured work directory: %s\n' "${#runs[@]}"
  printf 'Protected phase resume targets: %s\n' "${#protected[@]}"
}

cache_disk_usage() {
  printf 'category\tpath\tdisk_usage\n'
  local category path usage
  while IFS=$'\t' read -r category path; do
    if [[ -e "$path" ]]; then
      usage="$(du -sh "$path" 2>/dev/null | awk '{print $1}' || true)"
      [[ -n "$usage" ]] || usage="unavailable"
    else
      usage="0"
    fi
    printf '%s\t%s\t%s\n' "$category" "$path" "$usage"
  done <<EOF
nextflow_work	${WORK_DIR}
conda_cache	${CONDA_CACHE_DIR}
nextflow_metadata	${SCRIPT_DIR}/.nextflow/cache
publication_staging	${PUBLICATION_STAGING_DIR}
EOF
}

cache_prune_runs() {
  local runs=() protected=() i run_name retained_count=0 removed=0 failures=0 clean_output
  declare -A retained=()
  mapfile -t runs < <(cache_managed_runs)
  mapfile -t protected < <(cache_pointer_runs)
  for run_name in "${protected[@]:-}"; do
    [[ -n "$run_name" ]] || continue
    if cache_run_is_managed "$run_name" && nextflow log "$run_name" >/dev/null 2>&1; then
      retained["$run_name"]=1
    fi
  done
  if (( ${#retained[@]} == 0 )); then
    if adopt_legacy_cache_baseline; then
      run_name="$(head -n 1 "${RUNTIME_DIR}/last_successful_run.all" | tr -d '[:space:]')"
      retained["$run_name"]=1
    else
      echo "No protected phase resume target or valid legacy successful ASPIRE run exists; refusing to prune reusable work." >&2
      echo "Complete one run with the phase-aware controller before pruning." >&2
      exit 1
    fi
  fi
  if (( ${#retained[@]} > CACHE_KEEP )); then
    echo "Cannot retain ${CACHE_KEEP} runs because ${#retained[@]} phase resume targets are protected." >&2
    exit 1
  fi
  for ((i=${#runs[@]}-1; i>=0 && ${#retained[@]}<CACHE_KEEP; i--)); do
    retained["${runs[$i]}"]=1
  done
  for run_name in "${runs[@]}"; do
    [[ -z "${retained[$run_name]:-}" ]] || continue
    if [[ "$CACHE_FORCE" -eq 1 ]]; then
      if clean_output="$(nextflow clean "$run_name" -f -q 2>&1)"; then
        printf 'removed\t%s\n' "$run_name"
      else
        printf 'failed\t%s\t%s\n' "$run_name" "${clean_output//$'\n'/ }" >&2
        failures=$((failures + 1))
      fi
    else
      printf 'would_remove\t%s\n' "$run_name"
    fi
    removed=$((removed + 1))
  done
  retained_count="${#retained[@]}"
  if [[ "$CACHE_FORCE" -eq 1 ]]; then
    printf 'Processed %s removable run record(s); retained %s; failures %s.\n' "$removed" "$retained_count" "$failures"
    (( failures == 0 )) || return 1
  else
    printf 'Dry run: %s run record(s) would be removed and %s retained. Add --force to apply.\n' "$removed" "$retained_count"
  fi
}

cache_clear_runs() {
  local runs=() run_name pointer
  if [[ "$CACHE_CLEAR_ALL" -eq 1 ]]; then
    mapfile -t runs < <(cache_managed_runs)
  else
    if ! cache_run_is_managed "$CACHE_CLEAR_RUN"; then
      echo "Unknown Nextflow run for the configured work directory: ${CACHE_CLEAR_RUN}" >&2
      exit 1
    fi
    runs=("$CACHE_CLEAR_RUN")
  fi
  for run_name in "${runs[@]}"; do
    if [[ "$CACHE_FORCE" -eq 1 ]]; then
      nextflow clean "$run_name" -f -q
      printf 'removed\t%s\n' "$run_name"
      for pointer in "${RUNTIME_DIR}"/last_successful_run.*; do
        [[ -f "$pointer" ]] || continue
        [[ "$(head -n 1 "$pointer" | tr -d '[:space:]')" == "$run_name" ]] && rm -f "$pointer"
      done
    else
      printf 'would_remove\t%s\n' "$run_name"
    fi
  done
  if [[ "$CACHE_FORCE" -ne 1 ]]; then
    echo "Dry run only; add --force to remove the selected run cache/history."
  fi
  if [[ "$CACHE_CLEAR_ALL" -eq 1 && "$CACHE_INCLUDE_CONDA" -eq 1 ]]; then
    if [[ "$CACHE_FORCE" -ne 1 ]]; then
      echo "Conda cache would also be cleared."
    else
      conda_cache_real="$(realpath -m "$CONDA_CACHE_DIR")"
      if [[ "$conda_cache_real" == "/" || "$conda_cache_real" == "$HOME" || ${#conda_cache_real} -lt 12 ]]; then
        echo "Refusing unsafe Conda cache target: ${conda_cache_real}" >&2
        exit 1
      fi
      find "$conda_cache_real" -mindepth 1 -maxdepth 1 \
        ! -name '.aspire-run.lock' -exec rm -rf -- {} +
      echo "Cleared shared Conda environments and package cache: ${conda_cache_real}"
    fi
  fi
}

if [[ "$CACHE_MODE" -eq 1 ]]; then
  case "$CACHE_ACTION" in
    list|status) cache_list_runs ;;
    usage) cache_disk_usage ;;
    prune) cache_prune_runs ;;
    clear) cache_clear_runs ;;
  esac
  exit 0
fi

# One workflow may own a Conda cache at a time. This makes orphaned Nextflow
# environment markers safe to remove after an interrupted run.
stale_env_locks=()
while IFS= read -r -d '' stale_lock; do
  stale_env_locks+=("$stale_lock")
done < <(find "$CONDA_CACHE_DIR" -maxdepth 1 -type f -name '.env-*.lock' -print0)
if (( ${#stale_env_locks[@]} > 0 )); then
  rm -f "${stale_env_locks[@]}"
  echo "[controller] Removed ${#stale_env_locks[@]} stale Nextflow Conda environment lock marker(s)."
fi
if [[ "$RESUME_ENABLED" -eq 0 && "$ASPIRE_PHASE" != "analysis" && -d "$PUBLICATION_STAGING_DIR" ]]; then
  publication_staging_real="$(realpath -m "$PUBLICATION_STAGING_DIR")"
  runtime_dir_real="$(realpath -m "$RUNTIME_DIR")"
  if [[ "$publication_staging_real" == "$runtime_dir_real"/* ]]; then
    rm -rf "$publication_staging_real"
  else
    echo "Refusing to clear publication staging outside runtime directory: ${publication_staging_real}" >&2
    exit 1
  fi
fi
mkdir -p "$PUBLICATION_STAGING_DIR"
mkdir -p "${PUBLICATION_STAGING_DIR}/logs"
echo "[controller] Runtime directory: ${RUNTIME_DIR}"
echo "[controller] Publication staging directory: ${PUBLICATION_STAGING_DIR}"
echo "[controller] Execution phase: ${ASPIRE_PHASE}"

has_nextflow_arg() {
  local expected="$1"
  local arg
  for arg in "${NEXTFLOW_ARGS[@]}"; do
    if [[ "$arg" == "$expected" || "$arg" == "${expected}="* ]]; then
      return 0
    fi
  done
  return 1
}

DEFAULT_NEXTFLOW_REPORT_ARGS=()
if ! has_nextflow_arg -with-report; then
  DEFAULT_NEXTFLOW_REPORT_ARGS+=(-with-report "${PUBLICATION_STAGING_DIR}/logs/nextflow_report.html")
fi
if ! has_nextflow_arg -with-timeline; then
  DEFAULT_NEXTFLOW_REPORT_ARGS+=(-with-timeline "${PUBLICATION_STAGING_DIR}/logs/nextflow_timeline.html")
fi
if ! has_nextflow_arg -with-trace; then
  DEFAULT_NEXTFLOW_REPORT_ARGS+=(-with-trace "${PUBLICATION_STAGING_DIR}/logs/nextflow_trace.tsv")
fi
if ! has_nextflow_arg -with-dag; then
  DEFAULT_NEXTFLOW_REPORT_ARGS+=(-with-dag "${PUBLICATION_STAGING_DIR}/logs/nextflow_dag.html")
fi

# Avoid collisions when reusing staging after a failed or retained run.
rm -f \
  "${PUBLICATION_STAGING_DIR}/logs/nextflow_report.html" \
  "${PUBLICATION_STAGING_DIR}/logs/nextflow_timeline.html" \
  "${PUBLICATION_STAGING_DIR}/logs/nextflow_trace.tsv" \
  "${PUBLICATION_STAGING_DIR}/logs/nextflow_dag.html"

printf '%q ' "$SCRIPT_PATH" "${ORIGINAL_ARGS[@]}" \
  > "${PUBLICATION_STAGING_DIR}/logs/launch_command.txt"
printf '\n' >> "${PUBLICATION_STAGING_DIR}/logs/launch_command.txt"
nextflow -version > "${PUBLICATION_STAGING_DIR}/logs/nextflow_version.txt" 2>&1

export NXF_WORK="$WORK_DIR"
export NXF_CONDA_CACHEDIR="$CONDA_CACHE_DIR"
export CONDA_PKGS_DIRS="${CONDA_CACHE_DIR}/pkgs"
export CONDA_REMOTE_MAX_RETRIES="${CONDA_REMOTE_MAX_RETRIES:-5}"
export CONDA_REMOTE_BACKOFF_FACTOR="${CONDA_REMOTE_BACKOFF_FACTOR:-2}"
export CONDA_REMOTE_CONNECT_TIMEOUT_SECS="${CONDA_REMOTE_CONNECT_TIMEOUT_SECS:-20}"
export CONDA_REMOTE_READ_TIMEOUT_SECS="${CONDA_REMOTE_READ_TIMEOUT_SECS:-120}"

ASPIRE_REAL_MAMBA="$(command -v mamba)"
if [[ -z "$ASPIRE_REAL_MAMBA" || ! -x "$ASPIRE_REAL_MAMBA" ]]; then
  echo "Controller environment does not provide a usable mamba executable." >&2
  exit 1
fi
export ASPIRE_REAL_MAMBA
export ASPIRE_MAMBA_LOCK_FILE="${TMPDIR:-/tmp}/aspire-mamba-${UID}.lock"
export ASPIRE_MAMBA_LOCK_TIMEOUT="${ASPIRE_MAMBA_LOCK_TIMEOUT:-1800}"
export ASPIRE_MAMBA_BUILD_TIMEOUT="${ASPIRE_MAMBA_BUILD_TIMEOUT:-900}"
export PATH="${SCRIPT_DIR}/processes/controller/bin:${PATH}"

if [[ "$RESUME_ENABLED" -eq 1 ]]; then
  BASELINE_RUN="$(select_baseline_run)"
else
  BASELINE_RUN=""
fi
if [[ -n "$BASELINE_RUN" ]]; then
  echo "[controller] Baseline run for ${ASPIRE_PHASE} cache/history lookup: ${BASELINE_RUN}"
else
  echo "[controller] No prior Nextflow run history found."
fi

RERUN_FROM_CANONICAL=""
if [[ -n "$RERUN_FROM" ]]; then
  RERUN_FROM_UPPER="$(printf '%s' "$RERUN_FROM" | tr '[:lower:]' '[:upper:]')"
  RERUN_FROM_CANONICAL="${PROCESS_ALIASES[$RERUN_FROM_UPPER]:-$RERUN_FROM_UPPER}"
  start_idx=-1
  for i in "${!PROCESS_ORDER[@]}"; do
    if [[ "${PROCESS_ORDER[$i]}" == "$RERUN_FROM_CANONICAL" ]]; then
      start_idx=$i
      break
    fi
  done

  if [[ $start_idx -lt 0 ]]; then
    echo "Unknown process for --rerun-from: $RERUN_FROM" >&2
    echo "Use --list-stages to see valid process names." >&2
    exit 1
  fi

fi

# Persistent per-stage generations provide cacheable, scoped reruns. Advancing
# a generation changes only that stage's task hash; no cached work is deleted.
CACHE_GENERATION_FILE="${RUNTIME_DIR}/cache_generations.tsv"
declare -A CACHE_GENERATIONS=()
for process_name in "${PROCESS_ORDER[@]}"; do
  CACHE_GENERATIONS["$process_name"]=0
done
if [[ -f "$CACHE_GENERATION_FILE" ]]; then
  while IFS=$'\t' read -r process_name generation; do
    [[ -n "${CACHE_GENERATIONS[$process_name]+set}" ]] || continue
    [[ "$generation" =~ ^[0-9]+$ ]] || continue
    CACHE_GENERATIONS["$process_name"]="$generation"
  done < "$CACHE_GENERATION_FILE"
fi
if [[ -n "$RERUN_FROM_CANONICAL" ]]; then
  CACHE_GENERATIONS["$RERUN_FROM_CANONICAL"]=$((CACHE_GENERATIONS[$RERUN_FROM_CANONICAL] + 1))
fi

CACHE_GENERATION_TMP="${CACHE_GENERATION_FILE}.tmp.$$"
: > "$CACHE_GENERATION_TMP"
for process_name in "${PROCESS_ORDER[@]}"; do
  printf '%s\t%s\n' "$process_name" "${CACHE_GENERATIONS[$process_name]}" >> "$CACHE_GENERATION_TMP"
done
mv "$CACHE_GENERATION_TMP" "$CACHE_GENERATION_FILE"

RERUN_CONFIG="${RUNTIME_DIR}/cache_generations.config"
{
  echo 'process {'
  for process_name in "${PROCESS_ORDER[@]}"; do
    generation="${CACHE_GENERATIONS[$process_name]}"
    printf "  withName: '%s' { ext.aspire_cache_generation = %s }\n" "$process_name" "$generation"
  done
  echo '}'
} > "$RERUN_CONFIG"

# Keep process definitions physically modular while preserving ASPIRE's single
# resolved configuration binding. Nextflow evaluates the composed source from
# the repository root, so projectDir, relative paths, task hashing, and normal
# resume semantics are unchanged.
GENERATED_PIPELINE="${SCRIPT_DIR}/.asv_pipeline.generated.$$.nf"
cleanup_generated_pipeline() {
  rm -f "$GENERATED_PIPELINE"
}
trap cleanup_generated_pipeline EXIT
awk '1' \
  "${SCRIPT_DIR}/asv_pipeline.nf" \
  "${SCRIPT_DIR}/workflow/modules/preprocessing.nf" \
  "${SCRIPT_DIR}/workflow/modules/taxonomy_decontamination.nf" \
  "${SCRIPT_DIR}/workflow/modules/metadata_ecology.nf" \
  "${SCRIPT_DIR}/workflow/modules/indicators_diagnostics.nf" \
  "${SCRIPT_DIR}/workflow/modules/networks_reporting.nf" \
  > "$GENERATED_PIPELINE"
if [[ -n "$RERUN_FROM_CANONICAL" ]]; then
  echo "[controller] Rerun generation advanced only for ${RERUN_FROM_UPPER} (${RERUN_FROM_CANONICAL})"
  echo "[controller] Nextflow will rerun true dependency descendants; independent branches remain cacheable"
fi

if [[ "$RESUME_ENABLED" -eq 1 ]]; then
  if [[ -n "$BASELINE_RUN" ]]; then
    RESUME_ARGS=(-resume "$BASELINE_RUN")
  else
    RESUME_ARGS=(-resume)
  fi
else
  RESUME_ARGS=()
  echo "[controller] Resume disabled for this run (--no-resume)."
fi

format_nextflow_console() {
  # DSL2 reports fully qualified WORKFLOW:PROCESS names.  Shorten only the
  # displayed wrapper hierarchy; internal process identities remain unchanged
  # so resume hashes, reports, and --rerun-from names stay stable.
  sed -u \
    -e 's/ASPIRE_PREPROCESS:/PREP:/g' \
    -e 's/ASPIRE_ANALYSIS:/ANALYSIS:/g' \
    -e 's/RUN_METADATA_ANALYSES:/COMMUNITY:/g' \
    -e 's/RUN_GRAPH_NETWORK:/GRAPH:/g' \
    -e 's/RUN_ASV_MAG_NETWORK:/ASV_MAG:/g' \
    -e 's/RUN_MODULE_MAG_ANCHORS:/MAG_ANCHOR:/g' \
    -e 's/RUN_GROUP_GUILD_FUNCTION:/GUILD:/g'
}

echo "[controller] Console workflow aliases: PREP, ANALYSIS, COMMUNITY, GRAPH, ASV_MAG, MAG_ANCHOR, GUILD"

ASPIRE_PHASE="$ASPIRE_PHASE" \
ASPIRE_PIPELINE_CONFIG="$CONFIG_FILE" \
SPARK_PIPELINE_CONFIG="$CONFIG_FILE" \
nextflow run "$GENERATED_PIPELINE" \
  -c "$RERUN_CONFIG" \
  --params-file "$CONFIG_FILE" \
  --pipeline_config "$CONFIG_FILE" \
  "${RESUME_ARGS[@]}" \
  -w "$NXF_WORK" \
  "${DEFAULT_NEXTFLOW_REPORT_ARGS[@]}" \
  "${NEXTFLOW_ARGS[@]}" 2>&1 \
  | format_nextflow_console \
  | tee "${PUBLICATION_STAGING_DIR}/logs/controller.log"

if has_nextflow_arg -preview; then
  echo "[controller] Preview completed; resume pointers and published outputs were not modified."
  exit 0
fi

COMPLETED_RUN="$(nextflow log -q 2>/dev/null | tail -n 1 || true)"
if [[ -n "$COMPLETED_RUN" ]]; then
  phase_pointer_tmp="${RUNTIME_DIR}/last_successful_run.${ASPIRE_PHASE}.tmp.$$"
  printf '%s\n' "$COMPLETED_RUN" > "$phase_pointer_tmp"
  mv "$phase_pointer_tmp" "${RUNTIME_DIR}/last_successful_run.${ASPIRE_PHASE}"
  TASK_LOG_DIR="${PUBLICATION_STAGING_DIR}/logs/tasks"
  mkdir -p "$TASK_LOG_DIR"
  nextflow log "$COMPLETED_RUN" -f 'process,hash,workdir,status,exit,duration' \
    > "${PUBLICATION_STAGING_DIR}/logs/task_execution.tsv"
  while IFS=$'\t' read -r process_name task_hash task_workdir task_status task_exit task_duration; do
    [[ -n "$process_name" && -d "$task_workdir" ]] || continue
    safe_process="$(printf '%s' "$process_name" | sed 's/[^A-Za-z0-9._-]/_/g')"
    safe_hash="$(printf '%s' "$task_hash" | sed 's/[^A-Za-z0-9._-]/_/g')"
    destination="${TASK_LOG_DIR}/${safe_process}/${safe_hash}"
    mkdir -p "$destination"
    for log_name in .command.sh .command.out .command.err .command.log .command.trace .exitcode; do
      [[ -f "${task_workdir}/${log_name}" ]] || continue
      cp "${task_workdir}/${log_name}" "${destination}/${log_name#.}"
    done
  done < "${PUBLICATION_STAGING_DIR}/logs/task_execution.tsv"
fi

if [[ -f "${SCRIPT_DIR}/.nextflow.log" ]]; then
  cp "${SCRIPT_DIR}/.nextflow.log" "${PUBLICATION_STAGING_DIR}/logs/nextflow.log"
fi

python "${SCRIPT_DIR}/processes/output_layout/organize_outputs.py" \
  --staging-dir "$PUBLICATION_STAGING_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --config "$CONFIG_FILE"

if [[ "$PRUNE_CACHE_AFTER_RUN" -gt 0 ]]; then
  CACHE_KEEP="$PRUNE_CACHE_AFTER_RUN"
  CACHE_FORCE=1
  echo "[controller] Pruning completed run/cache records to ${CACHE_KEEP} retained runs."
  cache_prune_runs
fi

if [[ "${KEEP_RUNTIME_DIR,,}" == "true" ]]; then
  echo "[controller] Retaining Nextflow resume state: ${RUNTIME_DIR}"
else
  runtime_dir_real="$(realpath -m "$RUNTIME_DIR")"
  output_dir_real="$(realpath -m "$OUTPUT_DIR")"
  script_dir_real="$(realpath -m "$SCRIPT_DIR")"
  if [[ "$runtime_dir_real" == "/" ||
        "$output_dir_real" == "$runtime_dir_real" ||
        "$output_dir_real" == "$runtime_dir_real"/* ||
        "$script_dir_real" == "$runtime_dir_real" ||
        "$script_dir_real" == "$runtime_dir_real"/* ]]; then
    echo "Refusing to remove unsafe runtime directory: ${runtime_dir_real}" >&2
    exit 1
  fi
  rm -rf "$runtime_dir_real"
  echo "[controller] Removed runtime directory after successful publication: ${runtime_dir_real}"
fi
