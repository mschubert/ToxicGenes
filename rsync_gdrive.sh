#!/bin/sh

cd $(dirname $0)
PROJNAME=${PWD##*/}

rsync -auvr $@ --prune-empty-dirs \
    --exclude "candidates_mech/" \
    --exclude "data/" \
    --exclude "backup/" \
    --include "*/" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --exclude="*" \
    . ~/Documents/Results/"$PROJNAME"

cd -
