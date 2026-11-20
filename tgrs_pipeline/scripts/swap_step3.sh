#!/bin/bash
# Cross-run of Step 3 (decays): run <install> on an activation file produced by the OTHER build.
#   bash swap_step3.sh <new_run> <setup_run> <activation_run> <install_dir> [total_seconds]
#     new_run        name of the new run folder under $TGRS_ROOT/runs (must not exist)
#     setup_run      existing run of the SAME install as <install_dir>: source files, geometry and the
#                    photon / prompt-proton results are taken from it (copied or symlinked, never changed)
#     activation_run existing run whose CosmicProtons.activation.dat is used (the other build's)
#     install_dir    MEGAlib install used for the decays
#   DRY=1 prepares the folder and prints what would be submitted, without calling sbatch.
# Defaults: 50000 s (half of a full Step 3), 12 instances per node, 100 s per instance.
set -e
NEW=${1:?new run name}; SETUP=${2:?setup run}; ACT=${3:?activation run}; INSTALL=${4:?MEGAlib install dir}
TOTAL=${5:-${TOTAL:-50000}}
ROOT=${TGRS_ROOT:-/global/scratch/users/$USER/tgrs}
SCRIPTS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DEC_TIME=${DEC_TIME:-100}; DEC_CPUS=${DEC_CPUS:-12}; THROTTLE=${THROTTLE:-32}; WALL=${WALL:-03:00:00}
S=$ROOT/runs/$SETUP; D=$ROOT/runs/$NEW; AF=$ROOT/runs/$ACT/CosmicProtons.activation.dat
[ -d "$S" ]  || { echo "No setup run $S"; exit 1; }
[ -s "$AF" ] || { echo "No activation file $AF"; exit 1; }
[ -e "$D" ]  && { echo "$D already exists; choose a new name"; exit 1; }
mkdir -p "$D/logs"
# Everything at the top level of the setup run except the decay results, the activation file,
# step folders, logs and old outputs: small regular files are copied, symlinks and big files are linked.
for f in "$S"/*; do
  b=$(basename "$f")
  case "$b" in
    CosmicProtons.Decays.p1*|CosmicProtons.activation.dat|step3_*|logs|*.root|*.pdf) continue;;
  esac
  [ -d "$f" ] && [ ! -L "$f" ] && continue
  if [ -L "$f" ] || [ "$(stat -c %s "$f")" -gt 20000000 ]; then ln -s "$(readlink -f "$f")" "$D/$b"
  else cp -p "$f" "$D/$b"; fi
done
# The tasks copy $RUNDIR/inputs/. into their own folders, so the run needs its own inputs folder.
[ -d "$S/inputs" ] || { echo "No inputs folder in $S"; exit 1; }
cp -r "$S/inputs" "$D/inputs"
rm -f "$D/inputs/CosmicProtons.activation.dat"   # the tasks link the activation file from $RUNDIR
cp -p "$AF" "$D/CosmicProtons.activation.dat"
[ -s "$D/inputs/CosmicProtons.Step3.source" ] || { echo "inputs/CosmicProtons.Step3.source missing in $D"; exit 1; }
echo "Run folder $D prepared: setup from $SETUP, activation from $ACT, install $INSTALL"
echo "Source files that mention the activation file or absolute paths (check they point to $D):"
grep -nH -E 'ActivationFile|/global/' "$D"/*.source 2>/dev/null | head -20 || true
N3=$(( (TOTAL + DEC_CPUS*DEC_TIME - 1) / (DEC_CPUS*DEC_TIME) ))
echo "Step 3: $N3 tasks x $DEC_CPUS instances x ${DEC_TIME}s = $((N3*DEC_CPUS*DEC_TIME)) s"
export RUNDIR=$D MEGALIB_INSTALL=$INSTALL SCRIPTS_DIR DEC_TIME TARGET_TIME=$TOTAL
sb () { sbatch --parsable -o "$D/logs/%x_%A_%a.out" "$@"; }
if [ -n "$DRY" ]; then echo "DRY run: would submit step3 array 1-$N3 and merge"; exit 0; fi
s3=$(sb -J ${NEW}_step3 --array=1-$N3%$THROTTLE --cpus-per-task=$DEC_CPUS -t $WALL --export=ALL,STAGE=step3 "$SCRIPTS_DIR/savio_task.sh")
m3=$(sb -J ${NEW}_merge3 --dependency=afterany:$s3 --export=ALL,STAGE=step3 "$SCRIPTS_DIR/savio_merge.sh")
echo "step3=$s3 merge=$m3"
