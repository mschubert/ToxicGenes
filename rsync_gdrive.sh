#!/bin/sh

rsync -auvr \
    --include "orf" \
    --include "*.pdf" \
    --include "*.xlsx" \
    --exclude="*" \
    $(dirname $0) ~/Documents/Work/ORFdosage
