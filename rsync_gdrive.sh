#!/bin/sh

rsync -auvr --prune-empty-dirs \
    --exclude "data/" \
    --exclude "backup/" \
    --include "*/" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --exclude="*" \
    $(dirname $0) ~/Documents/Work/ORFdosage
