#!/bin/bash

# Builds the two libraries the develop-cosi-dc6 branch of MEGAlib needs but
# setup.sh does not provide: cfitsio, then HEALPix C++ (healpix_cxx).
# Both go into a shared fc_cosi directory, using the SAME system compiler the
# MEGAlib/ROOT/Geant4 build used (no module loads, so ROOT is not rebuilt).
#
# Submit via:
#   sbatch -o /global/scratch/users/$USER/healpix_deps-%j.out build_healpix_deps.sh

#SBATCH -J BuildHealpixDeps
#SBATCH --account=fc_cosi
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

set -e

DEPS=/global/home/groups/fc_cosi/megalib_deps_dc6
JOBS=${SLURM_CPUS_PER_TASK:-4}

umask 022
mkdir -p "$DEPS"
cd "$DEPS"

echo "Host: ${HOSTNAME}"
echo "Compiler: $(gcc --version | head -1)"
echo "Deps directory: $DEPS"

# ---------------------------------------------------------------- cfitsio
echo "==== cfitsio ===="
CFV=""
for v in 4.6.2 4.6.0 4.5.0 4.4.1; do
  if curl -fL -o cfitsio.tar.gz "https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-$v.tar.gz"; then
    CFV=$v; break
  fi
done
[ -n "$CFV" ] || { echo "ERROR: could not download any cfitsio tarball"; exit 1; }
echo "Using cfitsio $CFV"

rm -rf "cfitsio-$CFV" "$DEPS/cfitsio"
tar xzf cfitsio.tar.gz
cd "cfitsio-$CFV"
./configure --prefix="$DEPS/cfitsio" --disable-curl
make -j"$JOBS"
make install
cd "$DEPS"

export PKG_CONFIG_PATH="$DEPS/cfitsio/lib/pkgconfig:$PKG_CONFIG_PATH"
pkg-config --exists cfitsio || { echo "ERROR: cfitsio.pc not found after install"; exit 1; }
echo "cfitsio OK: $(pkg-config --modversion cfitsio)"

# ---------------------------------------------------------------- HEALPix C++
echo "==== HEALPix C++ ===="
curl -fL -o healpix.tar.gz "https://sourceforge.net/projects/healpix/files/Healpix_3.83/Healpix_3.83_2024Nov13.tar.gz/download" \
  || { echo "ERROR: could not download HEALPix"; exit 1; }

rm -rf Healpix_3.83
tar xzf healpix.tar.gz
cd Healpix_3.83

export HEALPIX="$(pwd)"
export FITSDIR="$DEPS/cfitsio/lib"
export FITSINC="$DEPS/cfitsio/include"
export CC=gcc CXX=g++

# Non-interactive configure of the C++ package (also builds libsharp), then build
./configure --auto=cxx
make -j"$JOBS" cpp-all

# ---------------------------------------------------------------- verify
export PKG_CONFIG_PATH="$HEALPIX/lib/pkgconfig:$DEPS/cfitsio/lib/pkgconfig:$PKG_CONFIG_PATH"
echo "==== Verification ===="
ls -l "$HEALPIX/include/healpix_base.h" || { echo "ERROR: healpix_base.h not installed"; exit 1; }
pkg-config --exists healpix_cxx || { echo "ERROR: healpix_cxx.pc not found"; exit 1; }
echo "healpix_cxx OK: $(pkg-config --modversion healpix_cxx)"
echo "cflags: $(pkg-config --cflags healpix_cxx)"

chmod -R g-w "$DEPS"

echo
echo "SUCCESS. For the MEGAlib rebuild use:"
echo "  export HEALPIX=$HEALPIX"
echo "  export PKG_CONFIG_PATH=$HEALPIX/lib/pkgconfig:$DEPS/cfitsio/lib/pkgconfig:\$PKG_CONFIG_PATH"
