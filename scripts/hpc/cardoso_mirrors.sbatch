#!/bin/bash
#SBATCH --job-name=cardoso-mirrors
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-00:30:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --mail-type=END,FAIL
# --mail-user is deliberately NOT a #SBATCH line: Slurm reads directives
# literally and does not expand env vars, so a `--mail-user=$VAR` directive
# silently mails nobody. submit_cardoso_mirrors.sh passes it on the CLI.
#
# Builds mirrors/profiles/<acc>.tsv and mirrors/profile_quantiles/<acc>.tsv
# from the Expression Atlas reprocessing of Cardoso-Moreira et al. 2019.
# Single core is correct: the work is one streaming pass over an 86 MB TSV
# plus 119 sorts of ~59k floats. Sized from that, not guessed.

set -euo pipefail
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_DIR"

# --- Site (overridable; no absolute path is committed, $WORK expands per user) ---
HPC_ENV_MODULE="${HPC_ENV_MODULE:-hpc-env/13.1}"
PYTHON_MODULE="${PYTHON_MODULE:-Python/3.11.3-GCCcore-13.1.0}"
CARDOSO_RAW="${CARDOSO_RAW:-${WORK:?WORK is unset; pass CARDOSO_RAW explicitly}/hvantk-datasets/raw/expression/cardoso_e-mtab-6814/raw}"

# --- Required ---
ACCESSION="${ACCESSION:?Set ACCESSION, e.g. E-MTAB-6814}"
REGISTRY="${REGISTRY:?Set REGISTRY to the genes.tsv shipped with this job}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to a writable output directory}"

CONFIG_XML="$CARDOSO_RAW/${ACCESSION}-configuration.xml"
MATRIX="$CARDOSO_RAW/${ACCESSION}-tpms.tsv"

# --- Prove the inputs and the output path before the expensive pass ---
for f in "$CONFIG_XML" "$MATRIX" "$REGISTRY" ; do
  [[ -r "$f" ]] || { echo "ERROR: not readable: $f" >&2; exit 1; }
done
mkdir -p "$OUT_DIR"
: > "$OUT_DIR/.writetest" 2>/dev/null || { echo "ERROR: $OUT_DIR not writable" >&2; exit 1; }
rm -f "$OUT_DIR/.writetest"

# --- Modules ---
if command -v ml >/dev/null 2>&1; then
  ml purge || true
  ml "$HPC_ENV_MODULE" || echo "WARNING: could not load $HPC_ENV_MODULE" >&2
  ml "$PYTHON_MODULE"  || echo "WARNING: could not load $PYTHON_MODULE" >&2
fi
# A cluster module must not shadow a project env's own packages.
unset PYTHONPATH

echo "== cardoso mirrors =="
echo "date:      $(date -Iseconds)"
echo "host:      $(hostname)"
echo "job:       ${SLURM_JOB_ID:-<none>}"
echo "python:    $(python3 -V 2>&1)"
echo "accession: $ACCESSION"
echo "matrix:    $MATRIX"
echo "out:       $OUT_DIR"
echo "matrix sha256: $(sha256sum "$MATRIX" | cut -d' ' -f1)"
echo "config sha256: $(sha256sum "$CONFIG_XML" | cut -d' ' -f1)"

python3 "$PROJECT_DIR/convert_cardoso_moreira.py" \
  --config "$CONFIG_XML" \
  --matrix "$MATRIX" \
  --registry "$REGISTRY" \
  --dataset "$ACCESSION" \
  --out-profiles  "$OUT_DIR/profiles_${ACCESSION}.tsv" \
  --out-quantiles "$OUT_DIR/profile_quantiles_${ACCESSION}.tsv"

echo "-- output --"
wc -l "$OUT_DIR/profiles_${ACCESSION}.tsv" "$OUT_DIR/profile_quantiles_${ACCESSION}.tsv"
sha256sum "$OUT_DIR/profiles_${ACCESSION}.tsv" "$OUT_DIR/profile_quantiles_${ACCESSION}.tsv"
echo "done: $(date -Iseconds)"
