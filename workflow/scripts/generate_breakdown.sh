#!/bin/bash

set -euo pipefail

nucflag breakdown \
-i <(zcat -f $(find /project/logsdon_shared/projects/HPRC/CenMAP_chrY/KO_working/nucflag_v1.0/final -name "*.bed" -not -name "*.renamed.bed") | grep chrY | grep -v -P "random|#") \
-t percent
