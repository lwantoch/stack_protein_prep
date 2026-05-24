#!/usr/bin/env bash
# =============================================================================
# submit_gaussian.sh  (Slurm Option B: job arrays)
#
# CESGA-adapted version:
#   - CPU count is controlled by Slurm (--cpus / --cpus-per-task)
#   - Do NOT forward -nproc / %NProcShared to Gaussian
#   - Worker script should also be CESGA-adapted
# =============================================================================

set -euo pipefail

# -----------------------------
# Defaults: discovery / array
# -----------------------------
RECURSIVE=0
PATTERN="*.com"
ROOT="."
MAXPAR=50

# -----------------------------
# Defaults: policy / filtering
# -----------------------------
OVERWRITE=0
RESTART=0
ONLY_MISSING=0
ONLY_FAILED=0

# -----------------------------
# Overrides forwarded to worker
# -----------------------------
METHOD=""
BASIS=""
DISP=""

# -----------------------------
# Slurm resource defaults
# -----------------------------
SLURM_PARTITION="compute"
SLURM_TIME="24:00:00"
SLURM_MEM="32G"
SLURM_CPUS=16

JOBNAME_SINGLE="gauss1"
JOBNAME_ARRAY="gaussA"

# -----------------------------
# Paths
# -----------------------------
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER="${script_dir}/run_gaussian.sh"

usage() {
  cat <<'EOF'
submit_gaussian.sh — submit Gaussian jobs to Slurm (single job or job array)

SINGLE MODE:
  ./submit_gaussian.sh path/to/job.com [options]

ARRAY MODE (recursive):
  ./submit_gaussian.sh -r [-p PATTERN] [ROOT] [options]

Discovery:
  -r                  Enable recursive search (array mode)
  -p PATTERN          find -name PATTERN (default: "*.com")
  ROOT                Root directory (default: .)

Policy / filtering:
  -o                  Overwrite policy: submit even if successful output exists
  --only-missing      Submit only if expected *.log is missing
  --only-failed       Submit only if *.log exists but NOT "Normal termination"
  --restart           Ask worker to restart from checkpoint if possible

Overrides forwarded to worker:
  --method M          Set/replace method/functional
  --basis B           Set/replace basis set
  -d DISP             Set/replace EmpiricalDispersion=DISP,
                      use "-d none" to remove EmpiricalDispersion

Array control:
  -P MAXPAR           Max parallel array tasks (default: 50)

Slurm resources:
  --partition NAME    default: compute
  --time HH:MM:SS     default: 24:00:00
  --mem 32G           default: 32G
  --cpus N            default: 16

Examples:
  ./submit_gaussian.sh calc/job.com --method B3LYP --basis def2SVP -d GD3BJ
  ./submit_gaussian.sh -r -p "*.com" step02_gaussian/ --only-failed --restart -P 80
EOF
}

abspath() {
  python3 - "$1" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
}

derive_log_path() {
  local com="$1"
  local dir base
  dir="$(dirname "$com")"
  base="$(basename "$com" .com)"
  echo "${dir}/${base}.log"
}

log_is_success() {
  local log="$1"
  [[ -f "$log" ]] && grep -q "Normal termination of Gaussian" "$log"
}

should_include() {
  local com="$1"
  local log
  log="$(derive_log_path "$com")"

  if [[ "$OVERWRITE" -eq 1 ]]; then
    return 0
  fi

  if [[ "$ONLY_MISSING" -eq 1 ]]; then
    [[ ! -f "$log" ]]
    return $?
  fi

  if [[ "$ONLY_FAILED" -eq 1 ]]; then
    [[ -f "$log" ]] || return 1
    log_is_success "$log" && return 1
    return 0
  fi

  log_is_success "$log" && return 1
  return 0
}

require_worker() {
  if [[ ! -f "$WORKER" ]]; then
    echo "ERROR: Worker script not found: $WORKER" >&2
    exit 2
  fi
  if [[ ! -x "$WORKER" ]]; then
    echo "ERROR: Worker script is not executable: $WORKER" >&2
    echo "Fix: chmod +x run_gaussian.sh" >&2
    exit 2
  fi
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -r) RECURSIVE=1; shift ;;
      -p) PATTERN="${2:?ERROR: -p needs a pattern}"; shift 2 ;;
      -P) MAXPAR="${2:?ERROR: -P needs an int}"; shift 2 ;;

      -o) OVERWRITE=1; shift ;;
      --restart) RESTART=1; shift ;;
      --only-missing) ONLY_MISSING=1; shift ;;
      --only-failed) ONLY_FAILED=1; shift ;;

      --method) METHOD="${2:?ERROR: --method needs a value}"; shift 2 ;;
      --basis) BASIS="${2:?ERROR: --basis needs a value}"; shift 2 ;;
      -d) DISP="${2:?ERROR: -d needs a value}"; shift 2 ;;

      --partition) SLURM_PARTITION="${2:?ERROR: --partition needs a value}"; shift 2 ;;
      --time) SLURM_TIME="${2:?ERROR: --time needs a value}"; shift 2 ;;
      --mem) SLURM_MEM="${2:?ERROR: --mem needs a value}"; shift 2 ;;
      --cpus) SLURM_CPUS="${2:?ERROR: --cpus needs an int}"; shift 2 ;;

      -h|--help) usage; exit 0 ;;
      --) shift; break ;;
      -*)
        echo "ERROR: unknown option: $1" >&2
        usage
        exit 2
        ;;
      *)
        ROOT="$1"
        shift
        ;;
    esac
  done
}

validate_args() {
  if [[ "$ONLY_MISSING" -eq 1 && "$ONLY_FAILED" -eq 1 ]]; then
    echo "ERROR: --only-missing and --only-failed are mutually exclusive" >&2
    exit 2
  fi

  if [[ ! "$SLURM_CPUS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --cpus must be integer, got: $SLURM_CPUS" >&2
    exit 2
  fi
  if [[ ! "$MAXPAR" =~ ^[0-9]+$ ]]; then
    echo "ERROR: -P must be integer, got: $MAXPAR" >&2
    exit 2
  fi

  require_worker
  mkdir -p "${script_dir}/logs"
}

build_exports_base() {
  local exports="ALL"
  exports+=",GAUSS_WORKER=${WORKER}"

  exports+=",GAUSS_OVERWRITE=${OVERWRITE}"
  exports+=",GAUSS_RESTART=${RESTART}"
  exports+=",GAUSS_ONLY_MISSING=${ONLY_MISSING}"
  exports+=",GAUSS_ONLY_FAILED=${ONLY_FAILED}"

  exports+=",GAUSS_METHOD=${METHOD}"
  exports+=",GAUSS_BASIS=${BASIS}"
  exports+=",GAUSS_DISP=${DISP}"

  echo "$exports"
}

submit_single() {
  local input="$1"
  local abs
  abs="$(abspath "$input")"

  if [[ ! -f "$abs" ]]; then
    echo "ERROR: input file not found: $abs" >&2
    exit 2
  fi

  if ! should_include "$abs"; then
    echo "Skip by policy (nothing to do): $abs"
    exit 0
  fi

  local exports
  exports="$(build_exports_base)"
  exports+=",GAUSS_INPUT=${abs}"

  sbatch --export="$exports" <<SBATCH_EOF
#!/usr/bin/env bash
#SBATCH --job-name=${JOBNAME_SINGLE}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH --output=${script_dir}/logs/%x_%j.out
#SBATCH --error=${script_dir}/logs/%x_%j.err
#SBATCH --requeue

set -euo pipefail
mkdir -p "${script_dir}/logs"

RUN_FLAGS=()
[[ "\${GAUSS_OVERWRITE:-0}" -eq 1 ]] && RUN_FLAGS+=("-o")
[[ "\${GAUSS_RESTART:-0}" -eq 1 ]] && RUN_FLAGS+=("--restart")
[[ "\${GAUSS_ONLY_MISSING:-0}" -eq 1 ]] && RUN_FLAGS+=("--only-missing")
[[ "\${GAUSS_ONLY_FAILED:-0}" -eq 1 ]] && RUN_FLAGS+=("--only-failed")

[[ -n "\${GAUSS_METHOD:-}" ]] && RUN_FLAGS+=("--method" "\${GAUSS_METHOD}")
[[ -n "\${GAUSS_BASIS:-}"  ]] && RUN_FLAGS+=("--basis"  "\${GAUSS_BASIS}")
[[ -n "\${GAUSS_DISP:-}"   ]] && RUN_FLAGS+=("-d"       "\${GAUSS_DISP}")

"\${GAUSS_WORKER}" "\${RUN_FLAGS[@]}" "\${GAUSS_INPUT}"
SBATCH_EOF
}

submit_array() {
  if [[ ! -d "$ROOT" ]]; then
    echo "ERROR: ROOT is not a directory: $ROOT" >&2
    exit 2
  fi

  local list
  list="${script_dir}/logs/gaussian_inputs_${USER}_$(date +%Y%m%d_%H%M%S).txt"
  mkdir -p "${script_dir}/logs"
  : > "$list"

  while IFS= read -r -d '' f; do
    local abs
    abs="$(abspath "$f")"
    if should_include "$abs"; then
      echo "$abs" >> "$list"
    fi
  done < <(find "$ROOT" -type f -name "$PATTERN" -print0 | sort -z)

  local n
  n="$(wc -l < "$list" | tr -d ' ')"
  if [[ "$n" -eq 0 ]]; then
    echo "Nothing to submit (after filtering). List: $list"
    exit 0
  fi

  echo "Submitting ARRAY job:"
  echo "  tasks        : $n"
  echo "  max parallel : $MAXPAR"
  echo "  list         : $list"

  local exports
  exports="$(build_exports_base)"
  exports+=",GAUSS_LIST=${list}"

  sbatch --array=1-"$n"%${MAXPAR} --export="$exports" <<SBATCH_EOF
#!/usr/bin/env bash
#SBATCH --job-name=${JOBNAME_ARRAY}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${SLURM_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH --output=${script_dir}/logs/%x_%A_%a.out
#SBATCH --error=${script_dir}/logs/%x_%A_%a.err
#SBATCH --requeue

set -euo pipefail
mkdir -p "${script_dir}/logs"

INPUT="\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "\${GAUSS_LIST}")"
[[ -n "\$INPUT" ]] || { echo "ERROR: empty input for array id \$SLURM_ARRAY_TASK_ID" >&2; exit 2; }

RUN_FLAGS=()
[[ "\${GAUSS_OVERWRITE:-0}" -eq 1 ]] && RUN_FLAGS+=("-o")
[[ "\${GAUSS_RESTART:-0}" -eq 1 ]] && RUN_FLAGS+=("--restart")
[[ "\${GAUSS_ONLY_MISSING:-0}" -eq 1 ]] && RUN_FLAGS+=("--only-missing")
[[ "\${GAUSS_ONLY_FAILED:-0}" -eq 1 ]] && RUN_FLAGS+=("--only-failed")

[[ -n "\${GAUSS_METHOD:-}" ]] && RUN_FLAGS+=("--method" "\${GAUSS_METHOD}")
[[ -n "\${GAUSS_BASIS:-}"  ]] && RUN_FLAGS+=("--basis"  "\${GAUSS_BASIS}")
[[ -n "\${GAUSS_DISP:-}"   ]] && RUN_FLAGS+=("-d"       "\${GAUSS_DISP}")

"\${GAUSS_WORKER}" "\${RUN_FLAGS[@]}" "\$INPUT"
SBATCH_EOF
}

parse_args "$@"
validate_args

if [[ "$RECURSIVE" -eq 1 ]]; then
  submit_array
else
  if [[ "$ROOT" == "." ]]; then
    echo "ERROR: no input file provided for single mode. Use -r for recursive mode." >&2
    usage
    exit 2
  fi
  submit_single "$ROOT"
fi
