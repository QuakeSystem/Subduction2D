#!/bin/bash

REMOTE="bert@login01:/scratch/tectonics/bert/Subduction2D/Subduction2D_SZU2019/Figures/subduction_movies/"
LOCAL="$HOME/SD/Subduction2D/Subduction2D_SZU2019/subduction_movies/"

mkdir -p "$LOCAL"

rsync -avP "$REMOTE" "$LOCAL"
EOF