#!/bin/bash
# Run the patched TGRS program on one build's simulation results.
# Run it from inside the run folder (it holds the *.p1.sim.gz files, TGRS.dat, TGRS.geo.setup):
#   cd tgrs/runs/old_smoke
#   MEGALIB_INSTALL=/global/home/groups/fc_cosi/MEGAlib_geant4.11.4 TGRS_BIN=<path to>/tgrs_bin_dc6/TGRS bash ../../scripts/run_tgrs.sh
# For the small smoke files this takes minutes on one core and can run on a login node. For the large
# runs submit it instead:  sbatch --export=ALL ../../scripts/run_tgrs.sh
# Output: TGRS_output.root (histograms TGRS, Sim, Difference, Component_0..2).
# NO_ACTIVATION=1 leaves out the decay (activation) component: use it when that component is empty
# or incomplete. The result is then TGRS_output_noact.root and any existing TGRS_output.root is kept
# as TGRS_output_withact.root.

#SBATCH -J TGRS
#SBATCH --account=fc_cosi
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

MEGALIB_INSTALL=${MEGALIB_INSTALL:?set MEGALIB_INSTALL}
TGRS_BIN=${TGRS_BIN:?set TGRS_BIN to the TGRS program built for this install}
source "$MEGALIB_INSTALL/bin/source-megalib.sh" || { echo "Could not source MEGAlib"; exit 1; }

NOACT=${NO_ACTIVATION:-0}
NEED="CosmicPhotons.p1.sim.gz CosmicProtonsPrompt.p1.sim.gz TGRS.dat TGRS.geo.setup"
[ "$NOACT" = 1 ] || NEED="$NEED CosmicProtons.Decays.p1.sim.gz"
for f in $NEED; do
  [ -s "$f" ] || { echo "ERROR: missing $f in $(pwd)"; exit 1; }
done
if [ "$NOACT" = 1 ] && [ -s TGRS_output.root ]; then mv -f TGRS_output.root TGRS_output_withact.root; fi
rm -f TGRS_output.root

BKG=(-b "Cosmic Photons"        CosmicPhotons.p1.sim.gz
     -b "Cosmic Protons Prompt" CosmicProtonsPrompt.p1.sim.gz)
[ "$NOACT" = 1 ] || BKG+=(-b "Cosmic Protons Activation" CosmicProtons.Decays.p1.sim.gz)

# no -p: the detector energy resolution and noise are applied
"$TGRS_BIN" -t TGRS.dat -g TGRS.geo.setup "${BKG[@]}"
status=$?
echo "TGRS exit code: $status"
[ -s TGRS_output.root ] || { echo "ERROR: no TGRS_output.root produced"; exit 1; }
OUT=TGRS_output.root
if [ "$NOACT" = 1 ]; then mv -f TGRS_output.root TGRS_output_noact.root; OUT=TGRS_output_noact.root; fi
ls -l "$OUT"
exit "$status"
