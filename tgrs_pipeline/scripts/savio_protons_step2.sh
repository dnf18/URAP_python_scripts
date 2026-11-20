#!/bin/bash
# Savio submit script for TGRS proton activation - Step 2 (Activation Buildup)
# Submit via: sbatch savio_protons_step2.sh
#
# Depends on Step 1 having already finished in the SAME --chdir directory:
# it needs every IsotopesLarge.dat.p1.incN.dat file Step 1 wrote there.
#
# Step 2 is a bookkeeping/activation calculation, not an event simulation --
# there is no "Run" line in CosmicProtons.Step2.source, just a single
# Activator block. It is NOT run with mcosima; one plain `cosima` process
# reads ALL of Step 1's isotope-production files and writes one
# CosmicProtons.activation.dat. Running it with multiple instances would
# just redo the same deterministic calculation N times, not speed it up.
#
# The stock CosmicProtons.Step2.source hardcodes only 5
# IsotopesLarge.dat.p1.incN.dat files. If Step 1 was run with a different
# thread count (e.g. -t 32), that file is stale/wrong, so this script
# regenerates it first via build_step2_source.py, which globs the actual
# directory for whatever IsotopesLarge.dat.p1.inc*.dat files Step 1 really
# produced and writes CosmicProtons.Step2.generated.source.
#
# How to submit (use the wrapper so you get prompted for your email):
#   cd into the directory that holds your Step 1 output, then run
#   bash submit_step.sh <this_script>.sh
# The job runs in whatever directory you submit from (no --chdir needed).

#SBATCH -J ProtonStep2
#SBATCH --account=fc_cosi
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --signal=2@60

echo "Starting Step 2: Activation buildup on host ${HOST}..."

MEGALIB_INSTALL=${MEGALIB_INSTALL:-/global/home/groups/fc_cosi/MEGAlib_geant4.11.4}
echo "MEGAlib install: $MEGALIB_INSTALL"
# Use the SAME install here that you used for Step 1.
source "$MEGALIB_INSTALL/bin/source-megalib.sh" || exit 1

# Path to your clone of the shared URAP_python_scripts repo (for
# build_step2_source.py). Override by exporting SCRIPTS_DIR before submitting.
SCRIPTS_DIR=${SCRIPTS_DIR:-/global/home/users/$USER/URAP_python_scripts}

echo "Regenerating CosmicProtons.Step2.generated.source from Step 1's actual output..."
python3 "$SCRIPTS_DIR/build_step2_source.py" CosmicProtons.Step2.source "$(pwd)"

# Remove any old result so a stale file cannot pass the check below.
rm -f CosmicProtons.activation.dat

# Cosima prints a very long per-volume isotope listing (about 10 GB in a test
# run), so keep only the last 300 lines in the slurm log.
cosima CosmicProtons.Step2.generated.source 2>&1 | tail -n 300
status=${PIPESTATUS[0]}
echo "Cosima exit code: $status"

if [ ! -s CosmicProtons.activation.dat ]; then
  echo "ERROR: Step 2 produced no CosmicProtons.activation.dat"
  exit 1
fi

echo "Step 2 complete. CosmicProtons.activation.dat is ready for Step 3."
exit "$status"
