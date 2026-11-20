#!/bin/bash
# Merge the array-task outputs of one stage into the run folder (symlinks + one stub file).
# Submitted by submit_chain.sh with STAGE=photons|step1|step3.
#SBATCH --account=fc_cosi
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
: "${RUNDIR:?}" "${STAGE:?}" "${SCRIPTS_DIR:?}"
python3 "$SCRIPTS_DIR/merge_task_outputs.py" "$RUNDIR" "$STAGE"
