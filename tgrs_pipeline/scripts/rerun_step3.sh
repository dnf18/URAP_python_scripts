#!/bin/bash
# Redo ONLY Step 3 (decays) + its merge for an existing run, in fresh folders.
#   bash rerun_step3.sh <run_name> <megalib_install_dir> [total_seconds]
# Defaults (from the pilot on dc6): 12 instances per node (~5.3 GB each), 100 s per instance.
# Earlier step3_NNNN folders are moved to step3_attempt<N>/ so nothing is lost.
set -e
RUN=${1:?run name}; INSTALL=${2:?MEGAlib install dir}; TOTAL=${3:-${TOTAL:-100000}}
ROOT=${TGRS_ROOT:-/global/scratch/users/$USER/tgrs}
SCRIPTS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DEC_TIME=${DEC_TIME:-100}; DEC_CPUS=${DEC_CPUS:-12}; THROTTLE=${THROTTLE:-32}; WALL=${WALL:-03:00:00}
D=$ROOT/runs/$RUN
[ -s "$D/CosmicProtons.activation.dat" ] || { echo "No activation.dat in $D"; exit 1; }
N3=$(( (TOTAL + DEC_CPUS*DEC_TIME - 1) / (DEC_CPUS*DEC_TIME) ))
cd "$D"
if ls -d step3_[0-9][0-9][0-9][0-9] >/dev/null 2>&1; then
  k=1; while [ -e "step3_attempt$k" ]; do k=$((k+1)); done
  mkdir "step3_attempt$k"; mv step3_[0-9][0-9][0-9][0-9] "step3_attempt$k"/
  echo "Moved earlier Step 3 folders to $D/step3_attempt$k"
fi
find . -maxdepth 1 -type l -name 'CosmicProtons.Decays.p1.inc*' -delete
rm -f CosmicProtons.Decays.p1.sim.gz
echo "Step 3: $N3 tasks x $DEC_CPUS instances x ${DEC_TIME}s = $((N3*DEC_CPUS*DEC_TIME)) s"
export RUNDIR=$D MEGALIB_INSTALL=$INSTALL SCRIPTS_DIR DEC_TIME TARGET_TIME=$TOTAL
sb () { sbatch --parsable -o "$D/logs/%x_%A_%a.out" "$@"; }
s3=$(sb -J ${RUN}_step3 --array=1-$N3%$THROTTLE --cpus-per-task=$DEC_CPUS -t $WALL --export=ALL,STAGE=step3 "$SCRIPTS_DIR/savio_task.sh")
m3=$(sb -J ${RUN}_merge3 --dependency=afterany:$s3 --export=ALL,STAGE=step3 "$SCRIPTS_DIR/savio_merge.sh")
echo "step3=$s3 merge=$m3"
