#!/usr/bin/env bash
# ==============================================================================
# run_gaussian.sh — runs a single Gaussian16 job (.com input) on CESGA
#
# CESGA-specific notes:
#   - Load module with: module load cesga/2020 g16/c1
#   - Do NOT use %NProcShared in Gaussian input (stripped automatically below)
#   - Do NOT use %mem in Gaussian input (stripped automatically below)
#   - Number of cores is controlled by sbatch/srun -c
# ==============================================================================

set -euo pipefail

# -------- options (defaults) --------
OVERWRITE=0
RESTART=0
ONLY_MISSING=0
ONLY_FAILED=0

# Optional patch overrides
METHOD=""
BASIS=""
DISP=""

# -------- CLI parsing --------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) OVERWRITE=1; shift ;;

    --restart) RESTART=1; shift ;;

    --only-missing) ONLY_MISSING=1; shift ;;

    --only-failed) ONLY_FAILED=1; shift ;;

    # On CESGA this should not be used
    -nproc)
      echo "ERROR: -nproc / %NProcShared must not be used on CESGA." >&2
      echo "       Control cores via sbatch -c or srun -c instead." >&2
      exit 2
      ;;

    --method) METHOD="${2:?ERROR: --method needs a value}"; shift 2 ;;

    --basis) BASIS="${2:?ERROR: --basis needs a value}"; shift 2 ;;

    -d) DISP="${2:?ERROR: -d needs a value (e.g. GD3BJ or none)}"; shift 2 ;;

    --) shift; break ;;

    -*) echo "ERROR: unknown option: $1" >&2; exit 2 ;;

    *) break ;;
  esac
done

INPUT="${1:?ERROR: missing input .com file}"
[[ -f "$INPUT" ]] || { echo "ERROR: input not found: $INPUT" >&2; exit 2; }

WORKDIR="$(cd "$(dirname "$INPUT")" && pwd)"
INPUT_ABS="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"
BASE="$(basename "$INPUT" .com)"
LOG="${WORKDIR}/${BASE}.log"

cd "$WORKDIR"

# -------- helper --------
log_success() { [[ -f "$LOG" ]] && grep -q "Normal termination of Gaussian" "$LOG"; }

# -------- filter/skip policy --------
if [[ "$ONLY_MISSING" -eq 1 && "$ONLY_FAILED" -eq 1 ]]; then
  echo "ERROR: --only-missing and --only-failed are mutually exclusive" >&2
  exit 2
fi

if [[ "$OVERWRITE" -eq 0 ]]; then
  if [[ "$ONLY_MISSING" -eq 1 ]]; then
    if [[ -f "$LOG" ]]; then
      echo "SKIP(--only-missing): $LOG exists"
      exit 0
    fi

  elif [[ "$ONLY_FAILED" -eq 1 ]]; then
    if [[ ! -f "$LOG" ]]; then
      echo "SKIP(--only-failed): $LOG missing"
      exit 0
    fi
    if log_success; then
      echo "SKIP(--only-failed): already successful"
      exit 0
    fi

  else
    if log_success; then
      echo "SKIP: already successful"
      exit 0
    fi
  fi
fi

# -------- modules / gaussian --------
module purge
module load cesga/2020 g16/c1
command -v g16 >/dev/null || { echo "ERROR: g16 not found in PATH" >&2; exit 2; }

# -------- scratch --------
# CESGA sets TMPDIR=/scratch/<jobid> per job; use that instead of /scratch/<user>
# which requires a pre-created directory. Fall back to /tmp only for local runs.
_SCRATCH_BASE="${TMPDIR:-/tmp}"
export GAUSS_SCRDIR="${_SCRATCH_BASE}/gauss_${SLURM_JOB_ID:-local}_$$"
mkdir -p "$GAUSS_SCRDIR"

INPUT_TO_RUN="$INPUT_ABS"
# Temp files created below (stripped/restart/patched) are cleaned up at exit.
trap 'rm -rf "$GAUSS_SCRDIR" "${WORKDIR}/${BASE}.stripped.com" "${WORKDIR}/${BASE}.gpu.com" "${WORKDIR}/${BASE}.patched.com" "${WORKDIR}/${BASE}.restart.com" 2>/dev/null; true' EXIT

# -------- strip CESGA-forbidden link0 directives (always) --------
# MCPB.py and antechamber may write %NProcShared / %Mem into the input.
# On CESGA, core count is controlled by Slurm (--cpus-per-task); both
# directives must be absent or they conflict with the scheduler.
TMP_STRIPPED="${WORKDIR}/${BASE}.stripped.com"
python3 - "$INPUT_TO_RUN" "$TMP_STRIPPED" <<'PY'
import sys
inp, outp = sys.argv[1], sys.argv[2]
out = []
for ln in open(inp, encoding="utf-8", errors="replace"):
    s = ln.lstrip()
    if s.lower().startswith("%nproc"):
        continue
    if s.lower().startswith("%mem="):
        continue
    if s.lower().startswith("%cpu="):
        continue
    if s.lower().startswith("%gpucpu="):
        continue
    out.append(ln)
open(outp, "w", encoding="utf-8").write("".join(out))
PY
INPUT_TO_RUN="$TMP_STRIPPED"

# -------- inject GPU link0 directives --------
# SLURM sets CUDA_VISIBLE_DEVICES to allocated GPU device IDs.
# Gaussian 16 requires %CPU and %GPUCPU when running on GPUs.
# This block is a no-op for CPU-only jobs (copies file unchanged).
_N_CPUS="${SLURM_CPUS_PER_TASK:-1}"
_CUDA_DEVS="${CUDA_VISIBLE_DEVICES:-}"
TMP_GPU="${WORKDIR}/${BASE}.gpu.com"
python3 - "$INPUT_TO_RUN" "$TMP_GPU" "$_N_CPUS" "$_CUDA_DEVS" <<'PY'
import sys, re, shutil

inp, outp, n_cpus_str, cuda_devs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
n_cpus = max(1, int(n_cpus_str))

gpu_devs = []
if cuda_devs and cuda_devs not in ("NoDevFiles", ""):
    gpu_devs = [x.strip() for x in cuda_devs.split(",") if x.strip().isdigit()]

if not gpu_devs:
    shutil.copy(inp, outp)
    sys.exit(0)

cpu_list = ",".join(str(i) for i in range(n_cpus))
gpu_list = ",".join(gpu_devs)

lines = open(inp, encoding="utf-8", errors="replace").read().splitlines(True)
out = []
injected = False

for ln in lines:
    if ln.strip() == "--Link1--":
        injected = False
    if not injected and re.match(r"^#", ln.lstrip()):
        j = len(out) - 1
        while j >= 0 and out[j].strip() != "--Link1--" and not out[j].lstrip().startswith("%"):
            j -= 1
        pos = j + 1
        out = out[:pos] + [f"%CPU={cpu_list}\n", f"%GPUCPU={gpu_list}={cpu_list}\n"] + out[pos:]
        injected = True
    out.append(ln)

open(outp, "w", encoding="utf-8").write("".join(out))
PY
INPUT_TO_RUN="$TMP_GPU"

CHK_GUESS="${WORKDIR}/${BASE}.chk"

# -------- stale-chk guard --------
# Gaussian behaviour: when %Chk=foo.chk references a pre-existing .chk file
# with a completed prior calculation, Gaussian silently prepends
# "Geom=AllCheck Guess=Read" to the route.  That mutated route is incompatible
# with input files that still carry a title/charge/atom block (the layout
# MCPB.py emits for large_mk and small_opt), and Gaussian then crashes inside
# Link 602 ("Wanted a floating point number… Found an end-of-line").
# A fresh, non-restart submission must therefore start from a clean chk.
# Only the chk whose basename matches the input is deleted; shared chks (e.g.
# small_fc reading small_opt.chk) are untouched because their basename differs.
if [[ "$RESTART" -ne 1 && -f "$CHK_GUESS" ]]; then
  echo "Removing stale chk before fresh run: $CHK_GUESS"
  rm -f "$CHK_GUESS"
fi

# -------- restart logic --------
if [[ "$RESTART" -eq 1 && -f "$CHK_GUESS" ]]; then
  TMP_RESTART="${WORKDIR}/${BASE}.restart.com"

  python3 - "$INPUT_TO_RUN" "$TMP_RESTART" <<'PY'
import sys, re

inp, outp = sys.argv[1], sys.argv[2]
lines = open(inp, "r", encoding="utf-8", errors="replace").read().splitlines(True)

i = 0
while i < len(lines) and (lines[i].lstrip().startswith("%") or lines[i].strip() == ""):
    i += 1

if i >= len(lines) or not lines[i].lstrip().startswith("#"):
    raise SystemExit("ERROR: route section (#...) not found")

route_start = i
route_lines = []
while i < len(lines) and lines[i].lstrip().startswith("#"):
    route_lines.append(lines[i].rstrip("\n"))
    i += 1

route_text = " ".join(r.strip() for r in route_lines)

if "Geom=AllCheck" not in route_text:
    route_text += " Geom=AllCheck"
if "Guess=Read" not in route_text:
    route_text += " Guess=Read"

route_text = re.sub(r"\s+", " ", route_text).strip()

out = []
out.extend(lines[:route_start])
out.append(route_text + "\n")

rest = lines[i:]
blank_blocks = 0
for ln in rest:
    out.append(ln)
    if ln.strip() == "":
        blank_blocks += 1
        if blank_blocks >= 2:
            break

open(outp, "w", encoding="utf-8").write("".join(out))
PY

  INPUT_TO_RUN="$TMP_RESTART"
fi

# -------- patch overrides --------
if [[ -n "${METHOD}${BASIS}${DISP}" ]]; then
  TMP_PATCH="${WORKDIR}/${BASE}.patched.com"

  python3 - "$INPUT_TO_RUN" "$TMP_PATCH" "$METHOD" "$BASIS" "$DISP" <<'PY'
import sys, re

inp, outp, method, basis, disp = sys.argv[1:6]
method, basis, disp = method.strip(), basis.strip(), disp.strip()

lines = open(inp, "r", encoding="utf-8", errors="replace").read().splitlines(True)

i = 0
while i < len(lines) and (lines[i].lstrip().startswith("%") or lines[i].strip() == ""):
    i += 1
if i >= len(lines) or not lines[i].lstrip().startswith("#"):
    raise SystemExit("ERROR: route section (#...) not found")

route_start = i
route_lines = []
while i < len(lines) and lines[i].lstrip().startswith("#"):
    route_lines.append(lines[i].rstrip("\n"))
    i += 1

route_text = " ".join(r.strip() for r in route_lines)
route_text = re.sub(r"\s+", " ", route_text).strip()

level_pat = re.compile(r"\b([A-Za-z][A-Za-z0-9\-\+]+)\s*/\s*([A-Za-z0-9][A-Za-z0-9\-\+\(\)\*]+)\b")

if method or basis:
    m = level_pat.search(route_text)
    if m:
        old_m, old_b = m.group(1), m.group(2)
        new_m = method or old_m
        new_b = basis or old_b
        route_text = level_pat.sub(f"{new_m}/{new_b}", route_text, count=1)
    else:
        if not (method and basis):
            raise SystemExit("ERROR: no METHOD/BASIS token found; to add you must provide BOTH --method and --basis")
        route_text = f"{route_text} {method}/{basis}".strip()

disp_pat = re.compile(r"\bEmpiricalDispersion\s*=\s*(\S+)\b", flags=re.IGNORECASE)
if disp:
    if disp.lower() == "none":
        route_text = disp_pat.sub("", route_text)
    else:
        if disp_pat.search(route_text):
            route_text = disp_pat.sub(f"EmpiricalDispersion={disp}", route_text, count=1)
        else:
            route_text = f"{route_text} EmpiricalDispersion={disp}".strip()

route_text = re.sub(r"\s+", " ", route_text).strip()

lines[route_start] = route_text + "\n"
del lines[route_start + 1:route_start + len(route_lines)]

open(outp, "w", encoding="utf-8").write("".join(lines))
PY

  INPUT_TO_RUN="$TMP_PATCH"
fi

# -------- status output --------
echo "Host:    $(hostname)"
echo "Workdir: $WORKDIR"
echo "Input:   $INPUT_TO_RUN"
echo "Log:     $LOG"
echo "Scratch: $GAUSS_SCRDIR"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-unset}"

# -------- run gaussian --------
g16 < "$INPUT_TO_RUN" > "$LOG"

# -------- post-run verification --------
if ! tail -n 80 "$LOG" | grep -q "Normal termination of Gaussian"; then
  echo "ERROR: Gaussian did not finish normally: $LOG" >&2
  exit 4
fi

echo "OK: Normal termination"

# -------- auto-formchk for MCPB checkpoint files --------
# MCPB.py -s 2 (Seminario) reads Cartesian Force Constants from small_opt.fchk.
# The small_opt job writes geometry to small_opt.chk; the subsequent small_fc
# frequency job appends force constants to the SAME small_opt.chk
# (%Chk=..._small_opt.chk in the fc input). We therefore re-run formchk after
# both jobs: after small_opt (geometry only) and again after small_fc (geometry
# + force constants). The second run overwrites the first, giving MCPB.py the
# complete fchk it needs.
if [[ "$BASE" == *_small_opt ]]; then
  CHK_FILE="${WORKDIR}/${BASE}.chk"
  FCHK_FILE="${WORKDIR}/${BASE}.fchk"
  if [[ -f "$CHK_FILE" ]]; then
    echo "Running formchk: $CHK_FILE → $FCHK_FILE"
    formchk "$CHK_FILE" "$FCHK_FILE"
    echo "OK: formchk complete"
  else
    echo "WARNING: $CHK_FILE not found — skipping formchk" >&2
  fi
fi

# small_fc writes force constants back into small_opt.chk — regenerate the
# fchk so it contains the Cartesian Force Constants MCPB.py -s 2 requires.
if [[ "$BASE" == *_small_fc ]]; then
  OPT_BASE="${BASE%_fc}_opt"
  CHK_FILE="${WORKDIR}/${OPT_BASE}.chk"
  FCHK_FILE="${WORKDIR}/${OPT_BASE}.fchk"
  if [[ -f "$CHK_FILE" ]]; then
    echo "Regenerating formchk after small_fc: $CHK_FILE → $FCHK_FILE"
    formchk "$CHK_FILE" "$FCHK_FILE"
    echo "OK: formchk complete (with force constants)"
  else
    echo "WARNING: $CHK_FILE not found — skipping post-fc formchk" >&2
  fi
fi
