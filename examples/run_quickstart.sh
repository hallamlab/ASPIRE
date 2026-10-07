#!/usr/bin/env bash
# Build the bundled installation fixture, run ASPIRE, and validate every benchmark check.
set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT="${PWD}/aspire-quickstart"
THREADS=4
while (($#)); do
  case "$1" in
    --output) OUTPUT="${2:?--output requires a directory}"; shift 2 ;;
    --threads) THREADS="${2:?--threads requires a positive integer}"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [--output DIRECTORY] [--threads N]"
      echo 'Creates dataset/, run.yml and results/ under DIRECTORY; resumes existing runs.'
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo '--threads must be a positive integer' >&2; exit 2; }
mkdir -p "$OUTPUT"
OUTPUT="$(cd "$OUTPUT" && pwd)"
if [[ ! -d "$OUTPUT/dataset" ]]; then
  python3 "$PROJECT_DIR/examples/mock_test/build_quickstart.py" --output "$OUTPUT/dataset"
fi
"$PROJECT_DIR/examples/configure_mock_run.sh" \
  --dataset "$OUTPUT/dataset" --output "$OUTPUT/results" \
  --config-out "$OUTPUT/run.yml" --threads "$THREADS"
"$PROJECT_DIR/run_asv_pipeline.sh" "$OUTPUT/run.yml"
"$PROJECT_DIR/examples/validate_mock_run.sh" \
  --dataset "$OUTPUT/dataset" --results "$OUTPUT/results" 2>&1 | tee "$OUTPUT/validation.txt"
