#!/bin/bash
# Usage: bash new_run.sh [-a] [-v <version_name>]
# Run from the project root (Subduction2D/)
#   -a          adjoint model: use the driver files AND Project/Manifest.toml from
#               Subduction2D_SZU2019_sens/ (default: Subduction2D_SZU2019/ + root Project/Manifest)
#   -v <name>   name of the run directory / job (default: timestamp, e.g. run_20261124_1530)
#   -h          show this help
#
# Everything the run needs (driver files, Project.toml, Manifest.toml) is copied into the run
# folder, and the job runs with --project=<run folder>, so a queued job is unaffected by later
# edits or Pkg.update()s.

# This script solves the issue that Eejit uses the input files as the model is started.
# This is an issue when there's a queue and older runs use the newly overwritten input files.

usage() { echo "Usage: bash new_run.sh [-a] [-v <version_name>]"; }

VERSION=""
ADJOINT=false
while getopts ":av:h" opt; do
    case $opt in
        a) ADJOINT=true ;;
        v) VERSION=$OPTARG ;;
        h) usage; exit 0 ;;
        :) echo "Error: -$OPTARG requires an argument"; usage; exit 1 ;;
        \?) echo "Error: unknown option -$OPTARG"; usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -gt 0 ]; then
    echo "Error: unexpected argument(s): $*  (did you mean -v $1 ?)"
    usage; exit 1
fi

VERSION=${VERSION:-run_$(date +%Y%m%d_%H%M)}

# Report which JustRelax is used, from the [[deps.JustRelax]] entry in the copied Manifest.toml.
# - local checkout (path = ...): its commit + uncommitted changes
# - GitHub (repo-url = ...): url, branch/tag, and git-tree-sha1 (fingerprint of the exact code)
justrelax_provenance() {
    local entry
    entry=$(awk '/^\[\[deps\.JustRelax\]\]/{f=1; next} /^\[\[/{f=0} f' "$1")
    local jr_path=$(echo "$entry" | sed -n 's/^path = "\(.*\)"/\1/p')
    if [ -n "$jr_path" ]; then
        echo "JustRelax source: path $jr_path"
        echo "JustRelax commit: $(git -C "$jr_path" rev-parse HEAD)"; git -C "$jr_path" status --porcelain
    else
        echo "JustRelax:"; echo "$entry" | grep -E '^(repo-url|repo-rev|git-tree-sha1|version) ='
    fi
}

if $ADJOINT; then
    TEMPLATE_DIR="Subduction2D_SZU2019_sens"   # driver files + Project/Manifest.toml
else
    TEMPLATE_DIR="Subduction2D_SZU2019"        # driver files; Project/Manifest.toml from root
fi
ENV_DIR=$($ADJOINT && echo "$TEMPLATE_DIR" || echo ".")
OUTDIR="$TEMPLATE_DIR/data/Subduction2D_JRv0.8.1_QS/$VERSION"
DRIVERS="Subduction2D.jl Subduction2D_setup.jl Subduction2D_rheology.jl"

for f in $DRIVERS; do
    [ -f "$TEMPLATE_DIR/$f" ] || { echo "Error: $TEMPLATE_DIR/$f not found"; exit 1; }
done
for f in Project.toml Manifest.toml; do
    [ -f "$ENV_DIR/$f" ] || { echo "Error: $ENV_DIR/$f not found (run Pkg.instantiate() there first)"; exit 1; }
done
if [ -e "$OUTDIR" ]; then
    echo "Error: $OUTDIR already exists, pick another -v name"; exit 1
fi

# Create run directory
mkdir -p "$OUTDIR"

# Copy driver files
for f in $DRIVERS; do cp "$TEMPLATE_DIR/$f" "$OUTDIR/"; done

# Copy the Julia environment. Relative `path = "..."` sources (e.g. ../JustRelax.jl) are made
# absolute, since they would otherwise be resolved relative to the run folder.
for f in Project.toml Manifest.toml; do
    ENV_ABS=$(cd "$ENV_DIR" && pwd) perl -pe \
      's{path = "(?!/)([^"]+)"}{require Cwd; "path = \"" . Cwd::abs_path("$ENV{ENV_ABS}/$1") . "\""}e' \
      "$ENV_DIR/$f" > "$OUTDIR/$f"
done

# Copy slurm template and fill in its placeholders: job name, the julia currently on your PATH
# (so it works for any user), and the run folder (copied driver + environment)
JULIA_BIN_DIR=$(dirname "$(command -v julia)")
sed \
  -e "s|__JOB_NAME__|$VERSION|g" \
  -e "s|__JULIA_BIN_DIR__|$JULIA_BIN_DIR|g" \
  -e "s|__RUN_DIR__|$OUTDIR|g" \
  run_eejit_sub2d.sbatch > "$OUTDIR/run_eejit_sub2d.sbatch"
if grep -q "__[A-Z_]*__" "$OUTDIR/run_eejit_sub2d.sbatch"; then
    echo "Error: unfilled placeholder(s) in $OUTDIR/run_eejit_sub2d.sbatch:"
    grep -n "__[A-Z_]*__" "$OUTDIR/run_eejit_sub2d.sbatch"; exit 1
fi

# Julia only reads LocalPreferences.toml next to the active Project.toml, so without this copy
# the run ignores the CUDA_Runtime_jll pin (12.9) and loads a different CUDA runtime.
# Use the environment's own preferences if it has them, else the repo root's.
if [ -f "$ENV_DIR/LocalPreferences.toml" ]; then
    cp "$ENV_DIR/LocalPreferences.toml" "$OUTDIR/"
elif [ -f LocalPreferences.toml ]; then
    cp LocalPreferences.toml "$OUTDIR/"
fi

cp -r LocalModules "$OUTDIR/"
{
  echo "date: $(date -Iseconds)"
  echo "user: $USER"
  echo "model: $($ADJOINT && echo adjoint || echo forward) ($TEMPLATE_DIR, environment from $ENV_DIR)"
  echo "Subduction2D commit: $(git rev-parse HEAD)"; git status --porcelain
  justrelax_provenance "$OUTDIR/Manifest.toml"
  echo "julia: $(command -v julia), $(julia --version)"
} > "$OUTDIR/provenance.txt"

echo "Output directory created: $OUTDIR"
echo "Input files are stored as:"
for f in $DRIVERS Project.toml Manifest.toml; do echo "  $OUTDIR/$f"; done
echo "Submitting job..."
sbatch "$OUTDIR/run_eejit_sub2d.sbatch"
