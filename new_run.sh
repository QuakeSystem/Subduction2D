#!/bin/bash
# Usage: bash new_run.sh <version_name>
# Run from the project root (Subduction2D/)

# This script solves the issue that Eejit uses the input files as the model is started.
# This is an issue when there's a queue and older runs use the newly overwritten input files.

VERSION=$1
if [ -z "$VERSION" ]; then
    echo "Usage: bash new_run.sh <version_name>"
    exit 1
fi

TEMPLATE_DIR="Subduction2D_SZU2019"
OUTDIR="Subduction2D_SZU2019/Figures/Subduction2D_nonuniform/$VERSION"

# Create run directory
mkdir -p "$OUTDIR"

# Copy templates
cp "$TEMPLATE_DIR/Subduction2D.jl"           "$OUTDIR/"
cp "$TEMPLATE_DIR/Subduction2D_setup.jl"     "$OUTDIR/"
cp "$TEMPLATE_DIR/Subduction2D_rheology.jl"             "$OUTDIR/"

# Copy slurm script and inject version name
sed \
  -e "s/--job-name .*/--job-name $VERSION/" \
  -e "s|julia --project \./Subduction2D_SZU2019/Subduction2D.jl|julia --project ./$OUTDIR/Subduction2D.jl|" \
  run_eejit_sub2d.sbatch > "$OUTDIR/run_eejit_sub2d.sbatch"

echo "Output directory created: $OUTDIR"
echo "Input files are stored as:"
echo "  $OUTDIR/Subduction2D.jl"
echo "  $OUTDIR/Subduction2D_setup.jl"
echo "  $OUTDIR/Subduction2D_rheology.jl"
echo "Submitting job..."
sbatch "$OUTDIR/run_eejit_sub2d.sbatch"