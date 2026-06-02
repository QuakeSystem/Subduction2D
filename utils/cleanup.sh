## Move all .png output files (zoom and full) to a subfolder.
# Run from a working directory containing multiple version output folders
# Working directory
#  |
# --- Working directory/v1
# --- Working directory/v2
# --- Working directory/v3

for dir in */; do
    cd "$dir"
    
    # full figures
    if ls fu*.png 2>/dev/null | grep -q .; then
        mkdir -p full && mv fu*.png full
    fi
    
    # zoom figures
    if ls zo*.png 2>/dev/null | grep -q .; then
        mkdir -p zoom && mv zo*.png zoom
    fi
    
    cd ..
done