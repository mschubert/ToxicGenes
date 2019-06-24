#!/bin/sh

rsync -auvr \
    --exclude ".snakemake" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --include "*/" \
    --exclude "backup" \
    --exclude="*" \
    rug:/data/p282396/ORFdosage/* $(dirname $0)
