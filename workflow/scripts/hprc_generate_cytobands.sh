#!/bin/bash

set -euo pipefail

bedtools subtract \
    -a <(grep chrY -h /project/logsdon_shared/projects/HPRC/CenMAP_chrY/results_*/2-concat_asm/*-asm-comb-dedup.fa.fai | awk -v OFS="\t" '{ print $1, 0, $2}') \
    -b /project/logsdon_shared/projects/HPRC/CenMAP_chrY/results_*/final/bed/*_complete_correct_cens.bed | \
awk -v OFS="\t" '{ print $1, $2, $3, "other", "gneg"}' | \
cat - \
    <(awk -v OFS="\t" '{
        len=$3-$2
        midpt=$2+int(len / 2)
        print $1, $2, midpt, "centromere", "acen"
        print $1, midpt, $3, "centromere", "acen"
    }' /project/logsdon_shared/projects/HPRC/CenMAP_chrY/results_*/final/bed/*_complete_correct_cens.bed) | \
sort -k1,1 -k2,2n -u
