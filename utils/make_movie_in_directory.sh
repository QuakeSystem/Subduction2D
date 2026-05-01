#!/bin/bash
# Make a movie of all .png files in the current working directory
module purge
module load tools
module load ffmpeg
FRAMERATE=${1:-10}
FOLDERNAME=$(basename "$PWD")
MOVIES_DIR="../movies"
mkdir -p "$MOVIES_DIR"
OUTPUT="$MOVIES_DIR/${FOLDERNAME}.mp4"

TMPDIR=$(mktemp -d)
MAX=$(ls *.png | grep -oP '\d+(?=\.png)' | sort -n | tail -1)
PADWIDTH=${#MAX}

for f in *.png; do
    NUM=$(echo "$f" | grep -oP '\d+(?=\.png)')
    PADDED=$(printf "%0${PADWIDTH}d" "$((10#$NUM))")
    PREFIX=$(echo "$f" | grep -oP '^[a-z_]+')
    ln -s "$(realpath "$f")" "$TMPDIR/${PREFIX}${PADDED}.png"
done

ffmpeg -framerate "$FRAMERATE" -pattern_type glob -i "$TMPDIR/*.png" -c:v libx264 -pix_fmt yuv420p "$OUTPUT"

rm -rf "$TMPDIR"
echo "Written: $OUTPUT"
