#!/bin/bash

# Builds MEGAlib branch "develop-cosi-dc6" (Geant4 11.4, full activation
# support) into the shared fc_cosi group directory, so the team can use it:
#   . /global/home/groups/fc_cosi/MEGAlib_geant4.11.4/bin/source-megalib.sh
#
# Multi-hour compile (Geant4 + ROOT + MEGAlib from source) -> sbatch job,
# not something to run on the login node.
#
# Submit via: sbatch build_megalib_dc6.sh
# Follow progress with: tail -f slurm-<jobid>.out   (run from the directory you submitted in)

#SBATCH -J BuildMEGAlibDC6
#SBATCH --account=fc_cosi
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --chdir=/global/home/groups/fc_cosi
#SBATCH -t 08:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --signal=2@60

set -e

echo "Starting build on host ${HOSTNAME}..."

BRANCH=develop-cosi-dc6
INSTALL_DIR=/global/home/groups/fc_cosi/MEGAlib_geant4.11.4

if [ -e "$INSTALL_DIR" ]; then
  echo "ERROR: $INSTALL_DIR already exists. Remove or rename it first so we don't build on top of a partial install."
  exit 1
fi

# umask 022 -> new files come out owner-writable only; the group can read/use.
umask 022
echo "umask is: $(umask)"

# Clone on the right branch AND pass --branch to setup.sh: otherwise setup.sh
# switches to the latest official release branch (Geant4 10.2) and you silently
# get the wrong version.
git clone --branch $BRANCH https://github.com/zoglauer/megalib.git "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo "Geant4 version this branch expects: $(cat config/AllowedGeant4Versions.txt)  (114 = Geant4 11.4)"

./setup.sh \
    --branch=$BRANCH \
    --maxthreads=$SLURM_CPUS_PER_TASK

echo "Build complete. Removing group-write (g-w) recursively..."
chmod -R g-w "$INSTALL_DIR"

echo "Done. Group members can now source:"
echo "  $INSTALL_DIR/bin/source-megalib.sh"
