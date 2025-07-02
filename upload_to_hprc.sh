
set -euo pipefail

SUBMISSION_ID="283198D5-BB8F-49CB-B1F9-CE1500812D8E"
SUBMISSION_NAME="HPRC_R2_CHRY_QC_NUCFLAG"
INPUT_DIR="/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/results/nucflag/final"

ssds staging upload \
--submission-id "${SUBMISSION_ID}" \
--name "${SUBMISSION_NAME}" \
"${INPUT_DIR}"
