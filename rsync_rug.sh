#!/bin/sh

rsync -auvr \
    --exclude "backup" \
    --include "orf" \
    --include "ccle" \
    --include "tcga" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --exclude="*" \
    rug:/data/p282396/ORFdosage ~/Documents/Work/ORFdosage
