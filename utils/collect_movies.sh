#!/bin/bash

SRC_ROOT="/scratch/tectonics/bert/Subduction2D/Subduction2D_SZU2019/Figures/Subduction2D_JRv0.5.1"
DST_ROOT="/scratch/tectonics/bert/Subduction2D/Subduction2D_SZU2019/Figures/subduction_movies"

mkdir -p "$DST_ROOT"

for version_dir in "$SRC_ROOT"/v*/; do
    version=$(basename "$version_dir")

    if [ -d "$version_dir/movies" ]; then
        mkdir -p "$DST_ROOT/$version"

        rsync -au \
            "$version_dir/movies/"*.mp4 \
            "$DST_ROOT/$version/" \
            2>/dev/null
    fi
done