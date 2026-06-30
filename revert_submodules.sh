#!/bin/bash
# Revert the submodule conversion applied by convert_to_submodules.py.
#
# Usage: bash revert_submodules.sh [dir1 dir2 ...]
#   Defaults to: src/ config_src/
#
# This script:
#   1. Deletes all *_s.F90 files created by the converter.
#   2. Restores the original module files from git.

set -euo pipefail

DIRS="${*:-src/ config_src/}"

# Find and delete submodule files
echo "Removing *_s.F90 submodule files..."
count=0
for dir in $DIRS; do
    while IFS= read -r -d '' f; do
        rm -f "$f"
        ((count++)) || true
    done < <(find "$dir" -name '*_s.F90' -print0)
done
echo "  Removed $count submodule file(s)."

# Restore original module files from git
echo "Restoring original module files from git..."
git checkout -- $DIRS
echo "Done."
