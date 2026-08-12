#!/bin/bash -l
#SBATCH --job-name=fms_unit_tests
#SBATCH --partition=rooster,chem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --mem=32GB
#SBATCH --output=unit_tests.%j.o
#SBATCH --error=unit_tests.%j.e
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load gcc/12.2

SRCDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
cd "$SRCDIR"

# Fast defaults. Override at submit time, e.g.:
#   sbatch --export=FORCE_CONFIGURE=1,UNITTEST_CLEAN=1 sub_unit_tests.sh
FORCE_CONFIGURE="${FORCE_CONFIGURE:-0}"
UNITTEST_CLEAN="${UNITTEST_CLEAN:-0}"
NPROC="${SLURM_CPUS_PER_TASK:-1}"

# Keep threaded libraries from oversubscribing cores during test execution.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

echo "=== Build started: $(date) ==="
echo "Host: $(hostname)"
echo "Dir:  $SRCDIR"
echo "CPUs: ${NPROC}"
echo "FORCE_CONFIGURE=${FORCE_CONFIGURE} UNITTEST_CLEAN=${UNITTEST_CLEAN}"
echo ""

# ── Step 1: configure (zero interface; no GPU/MPI needed for unit tests) ────
echo "--- configure.sh ---"
if [[ "${FORCE_CONFIGURE}" == "1" || ! -f "$SRCDIR/CONFIGFMS" ]]; then
    ./configure.sh zero
else
    echo "Skipping configure (CONFIGFMS present). Set FORCE_CONFIGURE=1 to force."
fi

# ── Step 2: build libfms.a ───────────────────────────────────────────────────
echo ""
echo "--- make makefmslib -j${NPROC} ---"
make makefmslib -j"${NPROC}"

# ── Step 3: build and run unit tests ─────────────────────────────────────────
echo ""
echo "--- unit tests ---"
cd "$SRCDIR/unit_tests"
if [[ "${UNITTEST_CLEAN}" == "1" ]]; then
    make clean
fi

FC_VAL="$(grep '^FC' "$SRCDIR/CONFIGFMS" | awk '{print $3}')"
FFLAGS_VAL="$(grep '^FFLAGS' "$SRCDIR/CONFIGFMS" | cut -d= -f2-)"

# Build tests in parallel, then run once.
make test -j"${NPROC}" FC="${FC_VAL}" FFLAGS="${FFLAGS_VAL}"
./test

echo ""
echo "=== Build finished: $(date) ==="
