
set -euo pipefail

SUBMISSION_ID="5E9D3123-C9D3-47D4-8BB7-D82FF0DC84EA"
SUBMISSION_NAME="HPRC_R2_QC_NUCFLAG"
INPUT_DIR="/project/logsdon_shared/projects/HPRC/CenMAP/nucflag_hprc/results_hprc/final"

ssds staging upload \
--submission-id "${SUBMISSION_ID}" \
--name "${SUBMISSION_NAME}" \
"${INPUT_DIR}"
