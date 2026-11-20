#!/bin/bash

# Stage 2: finish compiling MEGAlib (develop-cosi-dc6) now that HEALPix C++ and
# cfitsio are installed. This resumes the build that stopped at healpix_base.h:
# it does NOT rebuild ROOT or Geant4 and it does NOT re-clone MEGAlib.
# (setup.sh would have run exactly: source env file -> make -> move env file.)
#
# Submit via:
#   sbatch -o /global/scratch/users/$USER/finish_dc6-%j.out finish_megalib_dc6.sh

#SBATCH -J FinishMEGAlibDC6
#SBATCH --account=fc_cosi
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

INSTALL=/global/home/groups/fc_cosi/MEGAlib_geant4.11.4
DEPS=/global/home/groups/fc_cosi/megalib_deps_dc6
HEALPIXDIR=$DEPS/Healpix_3.83
CFITSIODIR=$DEPS/cfitsio
JOBS=${SLURM_CPUS_PER_TASK:-8}

echo "Host: ${HOSTNAME}"
echo "Compiler: $(gcc --version | head -1)"

# ---- HEALPix / cfitsio environment (what MEGAlib's Makefile.linuxgcc looks for)
export HEALPIX=$HEALPIXDIR
export PKG_CONFIG_PATH="$HEALPIXDIR/lib/pkgconfig:$CFITSIODIR/lib/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="$HEALPIXDIR/lib:$CFITSIODIR/lib:$LD_LIBRARY_PATH"

pkg-config --exists healpix_cxx || { echo "ERROR: pkg-config cannot find healpix_cxx"; exit 1; }
echo "healpix_cxx: $(pkg-config --modversion healpix_cxx)"
echo "cflags: $(pkg-config --cflags healpix_cxx)"
echo "libs:   $(pkg-config --libs healpix_cxx)"

# ---- MEGAlib environment written by setup.sh (the build stopped before it was
# ---- renamed to bin/source-megalib.sh)
if [ -f "$INSTALL/new-source-megalib.sh" ]; then
  ENVFILE="$INSTALL/new-source-megalib.sh"
elif [ -f "$INSTALL/bin/source-megalib.sh" ]; then
  ENVFILE="$INSTALL/bin/source-megalib.sh"
else
  echo "ERROR: no MEGAlib environment script found in $INSTALL"; exit 1
fi
echo "Using environment file: $ENVFILE"
source "$ENVFILE" || { echo "ERROR: could not source $ENVFILE"; exit 1; }

# env.sh may have reset these; set them again for the make
export HEALPIX=$HEALPIXDIR
export PKG_CONFIG_PATH="$HEALPIXDIR/lib/pkgconfig:$CFITSIODIR/lib/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="$HEALPIXDIR/lib:$CFITSIODIR/lib:$LD_LIBRARY_PATH"

cd "$MEGALIB" || { echo "ERROR: \$MEGALIB not set after sourcing"; exit 1; }
echo "Building in: $(pwd)"

make -j"$JOBS"
if [ "$?" != "0" ]; then
  echo "ERROR: make failed. Look above for the first 'error:' line."
  exit 1
fi

# ---- Finish what setup.sh would have done
if [ -f "$INSTALL/new-source-megalib.sh" ] && [ ! -f "$INSTALL/bin/source-megalib.sh" ]; then
  mv "$INSTALL/new-source-megalib.sh" "$INSTALL/bin/source-megalib.sh"
fi

# Make the runtime libraries findable whenever someone sources MEGAlib
if ! grep -q "megalib_deps_dc6" "$INSTALL/bin/source-megalib.sh"; then
  cat >> "$INSTALL/bin/source-megalib.sh" <<EOF

# HEALPix C++ and cfitsio (added after the dc6 build)
export HEALPIX=$HEALPIXDIR
export PKG_CONFIG_PATH="$HEALPIXDIR/lib/pkgconfig:$CFITSIODIR/lib/pkgconfig:\$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="$HEALPIXDIR/lib:$CFITSIODIR/lib:\$LD_LIBRARY_PATH"
EOF
fi

chmod -R g-w "$DEPS" "$INSTALL" 2>/dev/null

echo
echo "==== Verification ===="
for p in cosima mcosima revan mimrec; do
  if [ -x "$INSTALL/bin/$p" ]; then echo "OK      $p"; else echo "MISSING $p"; fi
done
echo "AllowedGeant4Versions: $(cat "$INSTALL/config/AllowedGeant4Versions.txt" 2>/dev/null)"

if [ -x "$INSTALL/bin/cosima" ]; then
  echo
  echo "SUCCESS. Use this in job scripts:"
  echo "  source $INSTALL/bin/source-megalib.sh"
else
  echo "ERROR: cosima was not built."
  exit 1
fi
