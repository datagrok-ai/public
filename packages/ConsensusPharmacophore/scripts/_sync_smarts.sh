#!/usr/bin/env bash
# Refresh files/pharmacophore-features.csv from packages/Chem/ (canonical source).
# Run from the package root: bash scripts/_sync_smarts.sh
set -euo pipefail
SRC="../Chem/files/pharmacophore-features.csv"
DST="files/pharmacophore-features.csv"
cp "$SRC" "$DST"
echo "Synced $DST from $SRC. Update files/pharmacophore-features.SOURCE.md with today's date."
