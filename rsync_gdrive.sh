#!/bin/sh

cd $(dirname $0)
PROJNAME=ToxicGenes
#PROJNAME=${PWD##*/}

rsync -auvr $@ --prune-empty-dirs \
    --exclude "candidates_mech/" \
    --exclude "data/" \
    --exclude "**/backup/" \
    --exclude "**/backup_tcga/" \
    --exclude "**/backup_ccle/" \
    --include "report/*.xlsx" \
    --exclude "**/*.xlsx" \
    --include "*/" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --exclude="*" \
    . ~/Documents/Projects/"$PROJNAME"/Analysis_output

cd -
