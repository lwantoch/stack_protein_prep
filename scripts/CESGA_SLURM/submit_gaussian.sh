#!/usr/bin/env bash
# =============================================================================
# submit_gaussian.sh
#
# CESGA Gaussian Slurm submitter.
#
# Important:
#   - --mem is total job memory.
#   - It does NOT multiply memory by CPU count.
#   - 4G with 32 CPUs means 4G total, not 125 GiB.
#   - CPU count is controlled by Slurm --cpus / --cpus-per-task.
#   - GPU count is controlled by --gpus.
#   - The actual Gaussian execution is delegated to run_gaussian.sh.
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
SLURM_PARTITION="short"
SLURM_TIME="06:00:00"
SLURM_MEM="4G"
SLURM_CPUS=32
SLURM_GPUS=1
SLURM_DEPENDENCY=""

JOBNAME_SINGLE="gauss1"
JOBNAME_ARRAY="gaussA"

# -----------------------------
# Paths
# -----------------------------
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER="${script_dir}/run_gaussian.sh"

usage() {
  cat <<'EOF'
submit_gaussian.sh — submit Gaussian jobs to Slurm

SINGLE MODE:
  ./submit_gaussian.sh path/to/job.com [options]

ARRAY MODE:
  ./submit_gaussian.sh -r [-p PATTERN] [ROOT] [options]

Discovery:
  -r                  Enable recursive search, submitting matching .com files as an array
  -p PATTERN          find -name PATTERN, default: "*.com"
  ROOT                Root directory for array mode, default: .

Policy / filtering:
  -o                  Submit even if successful output already exists
  --only-missing      Submit only if expected *.log is missing
  --only-failed       Submit only if *.log exists but does not contain normal termination
  --restart           Ask worker to restart from checkpoint if possible

Overrides forwarded to worker:
  --method M          Set/replace method/functional
  --basis B           Set/replace basis set
  -d DISP             Set/replace EmpiricalDispersion=DISP
                      use "-d none" to remove EmpiricalDispersion

Array control:
  -P MAXPAR           Max parallel array tasks, default: 50

Slurm resources:
  --partition NAME    default: short
  --time HH:MM:SS     default: 06:00:00
  --mem MEM           total job memory, default: 4G
  --cpus N            CPUs per task, default: 32
  --cpus-per-task N   alias for --cpus
  --gpus N            A100 GPUs, default: 1
  --dependency DEP    Slurm dependency, for example afterok:123456

Examples:
  ./submit_gaussian.sh calc/job.com
  ./submit_gaussian.sh calc/job.com --mem 4G --cpus 32 --gpus 1
  ./submit_gaussian.sh calc/job.com --dependency afterok:123456
  ./submit_gaussian.sh -r -p "*.com" step02_gaussian --only-failed --restart -P 20
EOF
}

abspath() {
  python3 - "$1" <<'PY'
import os
import sys

print(os.path.abspath(sys.argv[1]))
PY
}

derive_log_path() {
  local com="$1"
  local dir
  local base

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
    echo "ERROR: worker script not found: $WORKER" >&2
    exit 2
  fi

  if [[ ! -x "$WORKER" ]]; then
    echo "ERROR: worker script is not executable: $WORKER" >&2
    echo "Fix: chmod +x \"$WORKER\"" >&2
    exit 2
  fi
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -r)
        RECURSIVE=1
        shift
        ;;
      -p)
        PATTERN="${2:?ERROR: -p needs a pattern}"
        shift 2
        ;;
      -P)
        MAXPAR="${2:?ERROR: -P needs an int}"
        shift 2
        ;;

      -o)
        OVERWRITE=1
        shift
        ;;
      --restart)
        RESTART=1
        shift
        ;;
      --only-missing)
        ONLY_MISSING=1
        shift
        ;;
      --only-failed)
        ONLY_FAILED=1
        shift
        ;;

      --method)
        METHOD="${2:?ERROR: --method needs a value}"
        shift 2
        ;;
      --basis)
        BASIS="${2:?ERROR: --basis needs a value}"
        shift 2
        ;;
      -d)
        DISP="${2:?ERROR: -d needs a value}"
        shift 2
        ;;

      --partition)
        SLURM_PARTITION="${2:?ERROR: --partition needs a value}"
        shift 2
        ;;
      --time)
        SLURM_TIME="${2:?ERROR: --time needs a value}"
        shift 2
        ;;
      --mem)
        SLURM_MEM="${2:?ERROR: --mem needs a value}"
        shift 2
        ;;
      --cpus|--cpus-per-task)
        SLURM_CPUS="${2:?ERROR: $1 needs an int}"
        shift 2
        ;;
      --gpus)
        SLURM_GPUS="${2:?ERROR: --gpus needs an int}"
        shift 2
        ;;
      --dependency)
        SLURM_DEPENDENCY="${2:?ERROR: --dependency needs a value}"
        shift 2
        ;;

      -h|--help)
        usage
        exit 0
        ;;
      --)
        shift
        break
        ;;
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

  if [[ ! "$SLURM_GPUS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --gpus must be integer, got: $SLURM_GPUS" >&2
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
  local exports

  exports="ALL"
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

build_gpu_sbatch_line() {
  if [[ "$SLURM_GPUS" -gt 0 ]]; then
    echo "#SBATCH --gres=gpu:a100:${SLURM_GPUS}"
  fi
}

print_resource_summary() {
  echo "Slurm resources:"
  echo "  partition : $SLURM_PARTITION"
  echo "  time      : $SLURM_TIME"
  echo "  memory    : $SLURM_MEM total"
  echo "  cpus/task : $SLURM_CPUS"
  echo "  gpus/task : $SLURM_GPUS"
}

submit_single() {
  local input="$1"
  local abs
  local exports
  local dep_args=()
  local gpu_line

  abs="$(abspath "$input")"

  if [[ ! -f "$abs" ]]; then
    echo "ERROR: input file not found: $abs" >&2
    exit 2
  fi

  if ! should_include "$abs"; then
    echo "Skip by policy: $abs"
    exit 0
  fi

  exports="$(build_exports_base)"
  exports+=",GAUSS_INPUT=${abs}"

  if [[ -n "${SLURM_DEPENDENCY:-}" ]]; then
    dep_args=(--dependency="${SLURM_DEPENDENCY}")
  fi

  gpu_line="$(build_gpu_sbatch_line)"

  echo "Submitting single Gaussian job:"
  echo "  input     : $abs"
  print_resource_summary

  sbatch "${dep_args[@]}" --export="$exports" <<SBATCH_EOF
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
${gpu_line}

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
  local list
  local n
  local exports
  local gpu_line

  if [[ ! -d "$ROOT" ]]; then
    echo "ERROR: ROOT is not a directory: $ROOT" >&2
    exit 2
  fi

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

  n="$(wc -l < "$list" | tr -d ' ')"

  if [[ "$n" -eq 0 ]]; then
    echo "Nothing to submit after filtering."
    echo "List: $list"
    exit 0
  fi

  exports="$(build_exports_base)"
  exports+=",GAUSS_LIST=${list}"

  gpu_line="$(build_gpu_sbatch_line)"

  echo "Submitting Gaussian array job:"
  echo "  tasks        : $n"
  echo "  max parallel : $MAXPAR"
  echo "  list         : $list"
  print_resource_summary

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
${gpu_line}

set -euo pipefail
mkdir -p "${script_dir}/logs"

INPUT="\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "\${GAUSS_LIST}")"

if [[ -z "\$INPUT" ]]; then
  echo "ERROR: empty input for array id \${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

RUN_FLAGS=()
[[ "\${GAUSS_OVERWRITE:-0}" -eq 1 ]] && RUN_FLAGS+=("-o")
[[ "\${GAUSS_RESTART:-0}" -eq 1 ]] && RUN_FLAGS+=("--restart")
[[ "\${GAUSS_ONLY_MISSING:-0}" -eq 1 ]] && RUN_FLAGS+=("--only-missing")
[[ "\${GAUSS_ONLY_FAILED:-0}" -eq 1 ]] && RUN_FLAGS+=("--only-failed")

[[ -n "\${GAUSS_METHOD:-}" ]] && RUN_FLAGS+=("--method" "\${GAUSS_METHOD}")
[[ -n "\${GAUSS_BASIS:-}"  ]] && RUN_FLAGS+=("--basis" "\${GAUSS_BASIS}")
[[ -n "\${GAUSS_DISP:-}"   ]] && RUN_FLAGS+=("-d" "\${GAUSS_DISP}")

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
