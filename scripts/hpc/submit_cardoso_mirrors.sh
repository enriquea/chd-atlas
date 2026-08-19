#!/bin/bash
# Submit wrapper for cardoso_mirrors.sbatch.
#
# Sets --mail-user from $SLURM_EMAIL_NOTIFICATION (never a literal address) and
# forwards job vars through the environment. Override anything inline:
#   ACCESSION=E-MTAB-6814 OUT_DIR=/scratch/me/out bash submit_cardoso_mirrors.sh
set -euo pipefail
cd "$(dirname "$0")"

# Measured on this cluster with `sinfo`, not assumed: all_cpu.p has a 2h limit
# and idle capacity. Override CPU_PARTITION elsewhere.
CPU_PARTITION="${CPU_PARTITION:-all_cpu.p}"
ACCESSION="${ACCESSION:-E-MTAB-6814}"
REGISTRY="${REGISTRY:-$PWD/genes.tsv}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to a writable output directory}"

if [[ -z "${SLURM_EMAIL_NOTIFICATION:-}" && -r "$HOME/.bashrc" ]]; then
  SLURM_EMAIL_NOTIFICATION="$(
    sed -n "s/^[[:space:]]*\\(export[[:space:]]\\+\\)\\?SLURM_EMAIL_NOTIFICATION=[\"']\\?\\([^\"'#[:space:]]\\+\\).*/\\2/p" \
      "$HOME/.bashrc" | tail -n1
  )"
fi
export SLURM_EMAIL_NOTIFICATION="${SLURM_EMAIL_NOTIFICATION:-}"
if [[ -n "$SLURM_EMAIL_NOTIFICATION" ]]; then
  echo 'mail: ON  (address from $SLURM_EMAIL_NOTIFICATION; not printed)'
else
  echo 'mail: OFF (export SLURM_EMAIL_NOTIFICATION=you@example.org to enable)'
fi

export ACCESSION REGISTRY OUT_DIR
[[ -n "${CARDOSO_RAW:-}" ]] && export CARDOSO_RAW

SBATCH_ARGS=(--partition="$CPU_PARTITION" --export=ALL)
[[ -n "$SLURM_EMAIL_NOTIFICATION" ]] && SBATCH_ARGS+=(--mail-user="$SLURM_EMAIL_NOTIFICATION")

echo "partition: $CPU_PARTITION"
echo "accession: $ACCESSION"
echo "out dir:   $OUT_DIR"
# EXTRA_SBATCH_ARGS is optional. Note `"${ARR[@]:-}"` on an unset array expands
# to one EMPTY STRING argument, which sbatch then tries to open as a filename
# ("Unable to open file") -- so build a real array and expand it guardedly.
EXTRA=()
if [[ -n "${EXTRA_SBATCH_ARGS:-}" ]]; then
  # shellcheck disable=SC2206  # deliberate word-splitting: caller passes a flag string
  EXTRA=($EXTRA_SBATCH_ARGS)
fi

sbatch "${SBATCH_ARGS[@]}" ${EXTRA[@]+"${EXTRA[@]}"} cardoso_mirrors.sbatch
