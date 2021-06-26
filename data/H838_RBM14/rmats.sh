#!/usr/bin/env bash
python ~/rmats-turbo/rmats.py \
    --b1 luc.txt \
    --b2 rbm.txt \
    --gtf ref_annot.gtf \
    --readLength 100 \
    --nthread 12 \
    --od rmats_out \
    --tmp rmats_tmp

#-t paired \
