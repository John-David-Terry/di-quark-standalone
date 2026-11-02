#!/bin/zsh
# Download Pythia 8.311, apply the SRC patch (SRC_patch_v2.diff) and build it.
# Usage: setup_pythia.sh [install dir]
#   install dir defaults to $PYTHIA8_DIR, else ~/Data/pythia8311-SRCpatch (keep it outside the repo
#   and outside any cloud-synced folder). Set NCORES to limit the parallel build.
# The SRC:* settings are registered by the driver (fig3_patched.cc), not by this patch.
set -e
HERE=${0:A:h}
PYTHIA8_DIR=${1:-${PYTHIA8_DIR:-$HOME/Data/pythia8311-SRCpatch}}
URL=https://pythia.org/releases/pythia83/pythia8311.tgz
PATCH=$([[ -x /usr/bin/patch ]] && echo /usr/bin/patch || echo patch)
NCORES=${NCORES:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}

if [[ -e $PYTHIA8_DIR ]]; then
  echo "$PYTHIA8_DIR already exists; remove it or pass another directory" >&2
  exit 1
fi
TMP=$(mktemp -d)
trap 'rm -rf $TMP' EXIT

echo "downloading $URL"
curl -fL -o $TMP/pythia8311.tgz $URL
tar -xzf $TMP/pythia8311.tgz -C $TMP
mkdir -p ${PYTHIA8_DIR:h}
mv $TMP/pythia8311 $PYTHIA8_DIR
cd $PYTHIA8_DIR

# The patch is only valid on the pristine 8.311 versions of the three files it touches.
shasum -a 256 -c $HERE/pristine_8311.sha256
$PATCH -p1 < $HERE/SRC_patch_v2.diff
./configure > configure.log 2>&1
make -j$NCORES > build.log 2>&1
echo "built patched Pythia 8.311 in $PYTHIA8_DIR"
