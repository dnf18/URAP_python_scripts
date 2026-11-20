#!/bin/bash
# Build the patched MEGAlib TGRS program against ONE install.
#   bash build_tgrs.sh <megalib_install_dir> <that install's stock TGRS.cxx> <output_dir>
# e.g.
#   bash build_tgrs.sh /global/home/groups/fc_cosi/MEGAlib_geant4.11.4  inputs_dc6/TGRS.cxx tgrs_bin_dc6
#   bash build_tgrs.sh /global/home/groups/fc_cosi/MEGAlib_geant4.11.4 inputs_dc6/TGRS.cxx tgrs_bin_dc6
# Use each install's own TGRS.cxx and build it with that install: the two versions differ.
# Run it on a login node (one small compile), in a shell where no other MEGAlib is sourced.
INSTALL=${1:?megalib install dir}
SRC=${2:?stock TGRS.cxx of that install}
OUT=${3:?output dir}
HERE="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$OUT" || exit 1
python3 "$HERE/patch_tgrs.py" "$SRC" "$OUT/TGRS_patched.cxx" || exit 1

source "$INSTALL/bin/source-megalib.sh" || { echo "Could not source $INSTALL"; exit 1; }
LIBDIR="$MEGALIB/lib"
echo "MEGALIB=$MEGALIB"
echo "ROOT:   $(root-config --version)"

# MEGAlib ships either one combined library or several; link what is there.
if [ -f "$LIBDIR/libMEGAlib.so" ]; then
  LIBS="-lMEGAlib"
else
  LIBS=""
  for f in "$LIBDIR"/lib*.so; do b=$(basename "$f" .so); LIBS="$LIBS -l${b#lib}"; done
fi
echo "Linking: $LIBS"

g++ -O2 -o "$OUT/TGRS" "$OUT/TGRS_patched.cxx" \
    -I"$MEGALIB/include" $(root-config --cflags) \
    -L"$LIBDIR" -Wl,-rpath,"$LIBDIR" -Wl,--no-as-needed $LIBS $(root-config --glibs)
status=$?
if [ $status -ne 0 ]; then
  echo "BUILD FAILED. Send me the first 'error:' lines above."
  exit $status
fi
echo "Built $OUT/TGRS"
