#!/bin/sh

rsync -auvr --prune-empty-dirs \
    --include "*/" \
    --include "*.rds" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --exclude="*" \
    rug:/data/p282396/ORFdosage/* $(dirname $0)
