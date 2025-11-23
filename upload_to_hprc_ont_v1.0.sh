#!/bin/bash

set -euo pipefail

SUBMISSION_ID="283198D5-BB8F-49CB-B1F9-CE1500812D8E"
SUBMISSION_NAME="HPRC_R2_CHRY_QC_NUCFLAG"
INPUT_DIR="/project/logsdon_shared/projects/HPRC/CenMAP_chrY/KO_working/nucflag_v1.0/final"

wd=$(dirname $0)

ssds staging upload \
--submission-id "${SUBMISSION_ID}" \
--name "${SUBMISSION_NAME}" \
"${INPUT_DIR}" &> "${wd}/upload.log"

printf "sample_id\tdtype\turls\turi" > results_ont_v1.0.tsv
awk -v OFS="\t" -v FS="\t" '{
    match($1, "Copied (.+) to (.+)", arr);
    if (!arr[1]) { next };
    uri=arr[2]; url=arr[2];
    match(uri, "/([^_/]*?)/hprc_chry", sms);
    md5_suffix=match(uri, "md5") ? "_hash" : "";
    renamed_suffix=match(uri, "renamed") ? "_renamed" : "";
    file_desc="bedfile"renamed_suffix""md5_suffix
    gsub("s3://", "https://s3-us-west-2.amazonaws.com/", url);
    print sms[1], file_desc, url, uri
}' "${wd}/upload.log" | \
sort -k2,2 -k1,1 >> results_ont_v1.0.tsv
