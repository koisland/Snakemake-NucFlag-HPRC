#!/bin/bash

set -euo pipefail

nucflag ideogram \
    -i <(zcat -f $(find /project/logsdon_shared/projects/HPRC/CenMAP_chrY/KO_working/nucflag_v1.0/final -name "*.bed" -not -name "*.renamed.bed") | grep chrY) \
    -c /project/logsdon_shared/projects/Keith/NucFlag/exp/HPRC_chrY/cytobands.bed \
    -t 2 \
    -l 12_000_000
