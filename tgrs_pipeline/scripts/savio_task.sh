#!/bin/bash
# One array task of the TGRS chain. Do not sbatch this directly: submit_chain.sh does it.
# STAGE = photons | step1 | step3.  Each task runs in its own folder <RUNDIR>/<STAGE>_NNNN
# with a private copy of the inputs, so tasks never overwrite each other.
#SBATCH --account=fc_cosi
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1

: "${RUNDIR:?}" "${STAGE:?}" "${MEGALIB_INSTALL:?}"
echo "Host: $HOSTNAME  stage=$STAGE task=${SLURM_ARRAY_TASK_ID:-0}  install=$MEGALIB_INSTALL"
source "$MEGALIB_INSTALL/bin/source-megalib.sh" || { echo "Could not source MEGAlib"; exit 1; }
export OMP_NUM_THREADS=1
N=${SLURM_CPUS_PER_TASK:-1}

TD=$RUNDIR/${STAGE}_$(printf %04d "${SLURM_ARRAY_TASK_ID:-0}")
rm -rf "$TD"; mkdir -p "$TD"
cp -r "$RUNDIR/inputs/." "$TD/"
cd "$TD" || exit 1

setkey () {  # setkey <file> <key> <value> : replace the value of a key and check it
  sed -i -E "s/^(${2//./\\.}[[:space:]]+).*/\1${3}/" "$1"
  grep -E "^${2//./\\.}[[:space:]]+${3}\$" "$1" >/dev/null || { echo "ERROR: could not set $2 in $1"; exit 1; }
  echo "Set $2 = $3"
}

case "$STAGE" in
  photons)
    SRC=CosmicPhotons.source
    setkey $SRC SpaceSim.Time "${PHO_TIME:?}" ;;
  step1)
    SRC=CosmicProtons.Step1.source
    if grep -qE '^Level1\.Triggers' $SRC; then   # older MEGAlib stops on triggers
      setkey $SRC Level1.Triggers "${P1_TRIG:?}"
    else
      setkey $SRC Level1.Time "${P1_TIME:?}"
    fi ;;
  step3)
    SRC=CosmicProtons.Step3.source
    [ -s "$RUNDIR/CosmicProtons.activation.dat" ] || { echo "ERROR: no activation.dat from Step 2"; exit 1; }
    ln -s "$RUNDIR/CosmicProtons.activation.dat" .
    setkey $SRC ActivationStep3.Time "${DEC_TIME:?}" ;;
  *) echo "Unknown STAGE $STAGE"; exit 1 ;;
esac

mcosima -z -w -t "$N" -m "$N" "$SRC"
status=$?
echo "mcosima exit code: $status"
exit $status
