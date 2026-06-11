#!/usr/bin/env bash
set -euo pipefail

export NXF_VER="${NXF_VER:-25.10.0}"
export NXF_SYNTAX_PARSER="${NXF_SYNTAX_PARSER:-v1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "$0")"
CONTROL_ENV_DIR="${CONTROL_ENV_DIR:-${SCRIPT_DIR}/.controller_env}"

PROCESS_ORDER=(
  FASTP_QC
  MERGE_READS
  FILTER_READS
  RELABEL_FILTERED
  DEREPLICATE
  DENOISE
  CHIMERA_CHECK
  CREATE_COUNT_MATRIX
  FILTER_TABLE
  SINA_TRIM
  TAXONOMY
  MITOMASTER
  MITO_DECONTAM
  FILTER_COUNTS
  GENERAL_STATS
  PLOT_METADATA
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
  VOC_CORRELATION
  CLUSTERMAPS
  POWER_ANALYSIS_PIPELINE
  SPIECEASI
  NETWORK_MODULES
  ASV_MAG_LINK
  GRAPH_NETWORK
  MODULE_MAG_ANCHORS
  SANKEY
  MASTER_SUMMARY
)

usage() {
  cat <<'EOF'
Usage:
  run_asv_pipeline.sh [CONFIG_FILE] [--rerun-from PROCESS_NAME] [--resume-run RUN_NAME] [--resume-policy POLICY] [--no-resume] [--list-stages] [-- NEXTFLOW_ARGS...]

Options:
  --rerun-from PROCESS_NAME  Force rerun starting at PROCESS_NAME and all later processes in controller order, while preserving cacheability for future resumes.
  --resume-run RUN_NAME      Resume from a specific Nextflow run name/id (from `nextflow log -q`).
  --resume-policy POLICY     Resume baseline selection when --resume-run is not set.
                             Allowed: last-with-tasks (default), latest
  --no-resume                Disable resume for this run (cold execution).
  --list-stages              Print known process names for --rerun-from and exit.
  --help, -h                 Show this help.

Examples:
  run_asv_pipeline.sh asv_pipeline_nextflow.yml --rerun-from FILTER_COUNTS
  run_asv_pipeline.sh --resume-run lethal_poisson
  run_asv_pipeline.sh --resume-policy latest
  run_asv_pipeline.sh --rerun-from PLOT_METADATA -- -with-report report.html
EOF
}

if [[ -z "${IN_CONTROLLER_ENV:-}" ]]; then
  if ! command -v mamba >/dev/null 2>&1; then
    echo "mamba is required to bootstrap the controller environment." >&2
    exit 1
  fi
  if [[ ! -d "$CONTROL_ENV_DIR" ]]; then
    echo "[controller] Creating mamba env at $CONTROL_ENV_DIR"
    mamba env create --yes --prefix "$CONTROL_ENV_DIR" --file "${SCRIPT_DIR}/processes/controller/env.yml"
  fi
  exec env IN_CONTROLLER_ENV=1 CONTROL_ENV_DIR="$CONTROL_ENV_DIR" conda run --no-capture-output -p "$CONTROL_ENV_DIR" "$SCRIPT_PATH" "$@"
fi

CONFIG_FILE="asv_pipeline_nextflow.yml"
CONFIG_SET=0
LIST_STAGES=0
RERUN_FROM=""
NEXTFLOW_ARGS=()
RESUME_ENABLED=1
RESUME_POLICY="last-with-tasks"
RESUME_RUN=""
RERUN_CONFIG_FILE=""
declare -A RERUN_STAGE_SET=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --list-stages)
      LIST_STAGES=1
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
    --resume-policy)
      if [[ $# -lt 2 ]]; then
        echo "--resume-policy requires a value: last-with-tasks|latest" >&2
        exit 1
      fi
      RESUME_POLICY="$2"
      shift 2
      ;;
    --resume-policy=*)
      RESUME_POLICY="${1#*=}"
      shift
      ;;
    --no-resume)
      RESUME_ENABLED=0
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
      ;;
  esac
done

if [[ "$LIST_STAGES" -eq 1 ]]; then
  printf '%s\n' "${PROCESS_ORDER[@]}"
  exit 0
fi

if [[ "$RESUME_POLICY" != "last-with-tasks" && "$RESUME_POLICY" != "latest" ]]; then
  echo "Invalid --resume-policy '${RESUME_POLICY}'. Allowed: last-with-tasks, latest" >&2
  exit 1
fi

if [[ "$RESUME_ENABLED" -eq 0 && -n "$RESUME_RUN" ]]; then
  echo "--no-resume cannot be combined with --resume-run" >&2
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
        echo "[controller] Ignoring passthrough '-resume' argument; use --resume-run/--resume-policy/--no-resume." >&2
        ;;
      -resume=*)
        echo "[controller] Ignoring passthrough '-resume=...' argument; use --resume-run/--resume-policy/--no-resume." >&2
        ;;
      *)
        sanitized+=("$arg")
        ;;
    esac
  done
  NEXTFLOW_ARGS=("${sanitized[@]}")
}

has_tasks_for_run() {
  local run_name="$1"
  local rows rc process_name workdir _status
  set +e
  rows="$(nextflow log "$run_name" -f 'process,workdir,status' 2>&1)"
  rc=$?
  set -e
  [[ $rc -eq 0 ]] || return 1
  while IFS=$'\t' read -r process_name workdir _status; do
    [[ -n "$process_name" ]] || continue
    [[ "$process_name" == "process" && "$workdir" == "workdir" ]] && continue
    return 0
  done <<< "$rows"
  return 1
}

select_baseline_run() {
  local selected=""
  local latest=""
  local runs=()
  local i run_name

  if [[ -n "$RESUME_RUN" ]]; then
    if ! nextflow log "$RESUME_RUN" >/dev/null 2>&1; then
      echo "Specified --resume-run not found in Nextflow history: ${RESUME_RUN}" >&2
      exit 1
    fi
    selected="$RESUME_RUN"
    echo "$selected"
    return 0
  fi

  mapfile -t runs < <(nextflow log -q 2>/dev/null | sed '/^[[:space:]]*$/d')
  if [[ ${#runs[@]} -eq 0 ]]; then
    echo ""
    return 0
  fi
  latest="${runs[$(( ${#runs[@]} - 1 ))]}"

  if [[ "$RESUME_POLICY" == "latest" ]]; then
    echo "$latest"
    return 0
  fi

  for ((i=${#runs[@]}-1; i>=0; i--)); do
    run_name="${runs[$i]}"
    if has_tasks_for_run "$run_name"; then
      selected="$run_name"
      break
    fi
  done

  if [[ -n "$selected" ]]; then
    echo "$selected"
  else
    echo "$latest"
  fi
}

sanitize_nextflow_args

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  exit 1
fi

WORK_DIR=$(yq -r '.paths.work_dir // empty' "$CONFIG_FILE")
OUTPUT_DIR=$(yq -r '.paths.output_dir // empty' "$CONFIG_FILE")
CONDA_CACHE_DIR=$(yq -r '.paths.conda_cache_dir // empty' "$CONFIG_FILE")

if [[ -z "$WORK_DIR" ]]; then
  echo "paths.work_dir must be set in $CONFIG_FILE" >&2
  exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
  echo "paths.output_dir must be set in $CONFIG_FILE" >&2
  exit 1
fi

mkdir -p "$WORK_DIR"

if [[ -z "$CONDA_CACHE_DIR" || "$CONDA_CACHE_DIR" == "null" ]]; then
  CONDA_CACHE_DIR="${OUTPUT_DIR}/.conda_cache"
fi

mkdir -p "$CONDA_CACHE_DIR"

export NXF_WORK="$WORK_DIR"
export NXF_CONDA_CACHEDIR="$CONDA_CACHE_DIR"

BASELINE_RUN="$(select_baseline_run)"
if [[ -n "$BASELINE_RUN" ]]; then
  echo "[controller] Baseline run for cache/history lookup: ${BASELINE_RUN} (policy: ${RESUME_POLICY})"
else
  echo "[controller] No prior Nextflow run history found."
fi

if [[ -n "$RERUN_FROM" ]]; then
  RERUN_FROM_UPPER="$(printf '%s' "$RERUN_FROM" | tr '[:lower:]' '[:upper:]')"
  start_idx=-1
  for i in "${!PROCESS_ORDER[@]}"; do
    if [[ "${PROCESS_ORDER[$i]}" == "$RERUN_FROM_UPPER" ]]; then
      start_idx=$i
      break
    fi
  done

  if [[ $start_idx -lt 0 ]]; then
    echo "Unknown process for --rerun-from: $RERUN_FROM" >&2
    echo "Use --list-stages to see valid process names." >&2
    exit 1
  fi

  declare -A RERUN_STAGE_SET=()
  for ((i=start_idx; i<${#PROCESS_ORDER[@]}; i++)); do
    RERUN_STAGE_SET["${PROCESS_ORDER[$i]}"]=1
  done

  if [[ -z "$BASELINE_RUN" ]]; then
    echo "[controller] --rerun-from ${RERUN_FROM_UPPER}: no previous Nextflow run history found; proceeding without cache cleanup"
  else
    set +e
    task_rows="$(nextflow log "$BASELINE_RUN" -f 'process,workdir,status' 2>&1)"
    task_rows_rc=$?
    set -e
    if [[ $task_rows_rc -ne 0 ]]; then
      echo "[controller] Failed to inspect prior run tasks for cache cleanup (run: ${BASELINE_RUN})." >&2
      echo "[controller] nextflow log error:" >&2
      echo "${task_rows}" >&2
      echo "[controller] If another pipeline is currently running in this project, stop it and retry." >&2
      exit 1
    fi

    mapfile -t candidate_workdirs < <(
      while IFS=$'\t' read -r process_name workdir _status; do
        [[ -n "$process_name" ]] || continue
        process_key="$process_name"
        if [[ "$process_key" == *:* ]]; then
          process_key="${process_key##*:}"
        fi
        [[ -n "${RERUN_STAGE_SET[$process_name]:-}" || -n "${RERUN_STAGE_SET[$process_key]:-}" ]] || continue
        [[ -n "$workdir" && "$workdir" != "-" ]] || continue
        printf '%s\n' "$workdir"
      done <<< "$task_rows" | sort -u
    )

    if [[ ${#candidate_workdirs[@]} -eq 0 ]]; then
      echo "[controller] --rerun-from ${RERUN_FROM_UPPER}: no prior work directories found to invalidate"
    else
      work_root_real="$(realpath "$WORK_DIR")"
      removed_count=0
      for workdir in "${candidate_workdirs[@]}"; do
        workdir_real="$(realpath -m "$workdir")"
        if [[ "$workdir_real" == "$work_root_real"/* ]]; then
          rm -rf "$workdir_real"
          removed_count=$((removed_count + 1))
        else
          echo "[controller] Skipping unsafe work directory outside work root: ${workdir_real}" >&2
        fi
      done
      echo "[controller] Forcing rerun from ${RERUN_FROM_UPPER} by invalidating ${removed_count} prior task work directories"
      echo "[controller] Cache remains enabled; successful rerun tasks will be reusable on future -resume runs"
    fi
  fi
fi

if [[ -n "$RERUN_FROM" ]]; then
  RERUN_CONFIG_FILE="$(mktemp /tmp/aspire-rerun-cache-XXXXXX.config)"
  {
    echo "process {"
    for stage_name in "${PROCESS_ORDER[@]}"; do
      [[ -n "${RERUN_STAGE_SET[$stage_name]:-}" ]] || continue
      printf "  withName: /(^|.*:)%s\$/ { cache = false }\n" "$stage_name"
    done
    echo "}"
  } > "$RERUN_CONFIG_FILE"
  echo "[controller] Generated temporary Nextflow config to disable cache from ${RERUN_FROM_UPPER} onward: ${RERUN_CONFIG_FILE}"
fi

NEXTFLOW_RERUN_ARGS=()
if [[ -n "$RERUN_CONFIG_FILE" ]]; then
  NEXTFLOW_RERUN_ARGS=(-c "$RERUN_CONFIG_FILE")
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

ASPIRE_PIPELINE_CONFIG="$CONFIG_FILE" \
SPARK_PIPELINE_CONFIG="$CONFIG_FILE" \
nextflow run "${SCRIPT_DIR}/asv_pipeline.nf" \
  --params-file "$CONFIG_FILE" \
  --pipeline_config "$CONFIG_FILE" \
  "${RESUME_ARGS[@]}" \
  "${NEXTFLOW_RERUN_ARGS[@]}" \
  -w "$NXF_WORK" \
  "${NEXTFLOW_ARGS[@]}"
