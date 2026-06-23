#!/bin/bash

REMOTE="bert@eejit.geo.uu.nl:/scratch/tectonics/bert/Subduction2D/Subduction2D_SZU2019/data/subduction_movies/"
LOCAL="$HOME/SD/Subduction2D/Subduction2D_SZU2019/data/subduction_movies/"

mkdir -p "$LOCAL"

rsync -avP "$REMOTE" "$LOCAL"
EOF