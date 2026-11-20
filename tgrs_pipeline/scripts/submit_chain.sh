#!/bin/bash
# Submit the whole TGRS chain for one MEGAlib install as Slurm job arrays.
#   bash submit_chain.sh <run_name> <megalib_install_dir> <inputs_dir> [total_seconds]
# e.g. bash submit_chain.sh dc6_pilot /global/home/groups/fc_cosi/MEGAlib_geant4.11.4 ../inputs_dc6 5000
# total_seconds = simulated time per component (default 100000, the readme default).
# Needs, in the same folder as this script: savio_task.sh savio_merge.sh merge_task_outputs.py
#   savio_protons_step2.sh build_step2_source.py
# Chain: photons array | step1 array -> merge -> step2 -> step3 array -> merge
set -e
RUN=${1:?run name, e.g. dc6_pilot}
INSTALL=${2:?MEGAlib install dir}
INPUTS=$(cd "${3:?inputs dir (stock input files for this install)}" && pwd)
TOTAL=${4:-${TOTAL:-100000}}
ROOT=${TGRS_ROOT:-/global/scratch/users/$USER/tgrs}
SCRIPTS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Per-instance simulated seconds and CPUs per job (tuned to stay under the 3 h wall time).
P1_TIME=${P1_TIME:-18};   P1_CPUS=${P1_CPUS:-16}     # prompt protons (~7.7 s per instance-hour)
PHO_TIME=${PHO_TIME:-150}; PHO_CPUS=${PHO_CPUS:-32}  # photons (~70 s per instance-hour)
DEC_TIME=${DEC_TIME:-500}; DEC_CPUS=${DEC_CPUS:-25}  # decays (~600 s per instance-hour)
TRIG_RATE=${TRIG_RATE:-313.1}                         # old build: triggers per simulated second
THROTTLE=${THROTTLE:-32}                              # max array tasks running at once
WALL=${WALL:-03:00:00}

[ -f "$INSTALL/bin/source-megalib.sh" ] || { echo "No bin/source-megalib.sh in $INSTALL"; exit 1; }
for f in CosmicPhotons.source CosmicProtons.Step1.source CosmicProtons.Step2.source CosmicProtons.Step3.source; do
  [ -f "$INPUTS/$f" ] || { echo "Missing $f in $INPUTS"; exit 1; }
done
for f in savio_task.sh savio_merge.sh merge_task_outputs.py savio_protons_step2.sh build_step2_source.py; do
  [ -f "$SCRIPTS_DIR/$f" ] || { echo "Missing $f next to submit_chain.sh"; exit 1; }
done

ceil () { echo $(( ($1 + $2 - 1) / $2 )); }
N1=$(ceil "$TOTAL" $((P1_CPUS*P1_TIME)))
NP=$(ceil "$TOTAL" $((PHO_CPUS*PHO_TIME)))
N3=$(ceil "$TOTAL" $((DEC_CPUS*DEC_TIME)))
P1_TRIG=$(awk "BEGIN{printf \"%d\", $TRIG_RATE*$P1_TIME+0.5}")

D=$ROOT/runs/$RUN
[ -e "$D" ] && { echo "$D already exists; choose another run name"; exit 1; }
mkdir -p "$D/inputs" "$D/logs"
cp -r "$INPUTS"/. "$D/inputs/"
cp -r "$INPUTS"/. "$D/"            # Step 2 and TGRS read from the run folder itself
cd "$D"

export RUNDIR=$D MEGALIB_INSTALL=$INSTALL SCRIPTS_DIR P1_TIME P1_TRIG PHO_TIME DEC_TIME
echo "Allowed Geant4: $(cat "$INSTALL/config/AllowedGeant4Versions.txt" 2>/dev/null || echo unknown)"
echo "Install: $INSTALL   Inputs: $INPUTS   Run folder: $D"
echo "Target $TOTAL s per component:"
echo "  photons: $NP tasks x $PHO_CPUS instances x ${PHO_TIME}s"
echo "  step1  : $N1 tasks x $P1_CPUS instances x ${P1_TIME}s (old-build trigger cap $P1_TRIG if Level1.Triggers is used)"
echo "  step3  : $N3 tasks x $DEC_CPUS instances x ${DEC_TIME}s"

sb () { sbatch --parsable -o "$D/logs/%x_%A_%a.out" "$@"; }
pho=$(sb -J TGRSphotons --array=1-$NP%$THROTTLE --cpus-per-task=$PHO_CPUS -t $WALL --export=ALL,STAGE=photons "$SCRIPTS_DIR/savio_task.sh")
phm=$(sb -J merge_photons --dependency=afterany:$pho --export=ALL,STAGE=photons "$SCRIPTS_DIR/savio_merge.sh")
s1=$(sb -J TGRSstep1 --array=1-$N1%$THROTTLE --cpus-per-task=$P1_CPUS -t $WALL --export=ALL,STAGE=step1 "$SCRIPTS_DIR/savio_task.sh")
m1=$(sb -J merge_step1 --dependency=afterany:$s1 --export=ALL,STAGE=step1 "$SCRIPTS_DIR/savio_merge.sh")
s2=$(sb -J TGRSstep2 --dependency=afterok:$m1 -t 08:00:00 "$SCRIPTS_DIR/savio_protons_step2.sh")
s3=$(sb -J TGRSstep3 --dependency=afterok:$s2 --array=1-$N3%$THROTTLE --cpus-per-task=$DEC_CPUS -t $WALL --export=ALL,STAGE=step3 "$SCRIPTS_DIR/savio_task.sh")
m3=$(sb -J merge_step3 --dependency=afterany:$s3 --export=ALL,STAGE=step3 "$SCRIPTS_DIR/savio_merge.sh")
echo "photons=$pho merge=$phm | step1=$s1 merge=$m1 step2=$s2 step3=$s3 merge=$m3"
echo "Done when merge_photons and merge_step3 have finished: squeue -u $USER ; logs in $D/logs"
