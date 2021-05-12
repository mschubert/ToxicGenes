#!/usr/bin/env bash
python ~/rmats-turbo/rmats.py \
    --b1 luc.txt \
    --b2 rbm.txt \
    --gtf ./ref_annot.gtf \
    -t paired \
    --readLength 50 \
    --nthread 12 \
    --od rmats_out \
    --tmp rmats_tmp
