#!/bin/bash
# Re-submit the merge and the later steps of an existing run (when the tasks already finished
# but a merge failed). Does not rerun any simulation.
#   bash resume_chain.sh <run_name> <megalib_install_dir> [total_seconds]
# Before using it: scancel the old waiting jobs of that run (squeue -u $USER).
set -e
RUN=${1:?run name}
INSTALL=${2:?MEGAlib install dir}
TOTAL=${3:-${TOTAL:-100000}}
ROOT=${TGRS_ROOT:-/global/scratch/users/$USER/tgrs}
SCRIPTS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DEC_TIME=${DEC_TIME:-500}; DEC_CPUS=${DEC_CPUS:-25}; THROTTLE=${THROTTLE:-32}; WALL=${WALL:-03:00:00}
D=$ROOT/runs/$RUN
[ -d "$D" ] || { echo "No run folder $D"; exit 1; }
N3=$(( (TOTAL + DEC_CPUS*DEC_TIME - 1) / (DEC_CPUS*DEC_TIME) ))
cd "$D"
export RUNDIR=$D MEGALIB_INSTALL=$INSTALL SCRIPTS_DIR DEC_TIME
sb () { sbatch --parsable -o "$D/logs/%x_%A_%a.out" "$@"; }
mp=$(sb -J merge_photons --export=ALL,STAGE=photons "$SCRIPTS_DIR/savio_merge.sh")
m1=$(sb -J merge_step1 --export=ALL,STAGE=step1 "$SCRIPTS_DIR/savio_merge.sh")
s2=$(sb -J TGRSstep2 --dependency=afterok:$m1 -t 08:00:00 "$SCRIPTS_DIR/savio_protons_step2.sh")
s3=$(sb -J TGRSstep3 --dependency=afterok:$s2 --array=1-$N3%$THROTTLE --cpus-per-task=$DEC_CPUS -t $WALL --export=ALL,STAGE=step3 "$SCRIPTS_DIR/savio_task.sh")
m3=$(sb -J merge_step3 --dependency=afterany:$s3 --export=ALL,STAGE=step3 "$SCRIPTS_DIR/savio_merge.sh")
echo "merge_photons=$mp merge_step1=$m1 step2=$s2 step3=$s3 merge_step3=$m3"
