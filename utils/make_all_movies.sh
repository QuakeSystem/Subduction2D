#!/bin/bash
# Run from model_name directory
# Usage: bash ../../../../utils/make_all_movies.sh <framerate>
# Calls make_movie_in_directory.sh
FRAMERATE=${1:-10}
SCRIPT_DIR="$(dirname "$(realpath "$0")")"

for version_dir in */; do
    # Skip top-level movies folder
    [ "$version_dir" = "movies/" ] && continue

    for fig_dir in "${version_dir}fig_full" "${version_dir}fig_zoom" "${version_dir}full" "${version_dir}zoom"; do
        if [ -d "$fig_dir" ] && ls "$fig_dir"/*.png 2>/dev/null | grep -q .; then
            FOLDERNAME=$(basename "$fig_dir")
            MOVIE="${version_dir}movies/${FOLDERNAME}.mp4"
            if [ -f "$MOVIE" ]; then
                echo "Skipping $fig_dir (movie already exists)"
                continue
            fi
            echo "Making movie for: $fig_dir"
            cd "$fig_dir"
            bash "$SCRIPT_DIR/make_movie_in_directory.sh" "$FRAMERATE"
            cd - > /dev/null
        fi
    done
done