#!/bin/sh

cd $(dirname $0)
PROJNAME=${PWD##*/}

rsync -auvr $@ --prune-empty-dirs \
    --include "*/" \
    --include "*.rds" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --exclude="*" \
    rug:/data/p282396/"$PROJNAME"/* .

cd -
