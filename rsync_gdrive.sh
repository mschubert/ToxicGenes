#!/bin/sh

rsync -auvr \
    --exclude "backup" \
    --include "orf" \
    --include "ccle" \
    --include "tcga" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --exclude="*" \
    $(dirname $0) ~/Documents/Work/ORFdosage
