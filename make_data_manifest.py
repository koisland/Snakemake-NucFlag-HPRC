import os
import glob
import sys
import polars as pl

def main():
    INPUT_DIR = glob.glob("/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/results/*/final/*")
    SUBMISSION_ID="283198D5-BB8F-49CB-B1F9-CE1500812D8E"
    SUBMISSION_NAME="HPRC_R2_CHRY_QC_NUCFLAG"
    OUTPUT_URL = f"https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/{SUBMISSION_ID}--{SUBMISSION_NAME}/{{sample}}/hprc_chry/assembly_qc/nucflag/v0.3.4_hifi/"
    OUTPUT_URI = f"s3://human-pangenomics/submissions/{SUBMISSION_ID}--{SUBMISSION_NAME}/{{sample}}/hprc_chry/assembly_qc/nucflag/v0.3.4_hifi/"

    samples = []
    output_urls = []
    output_uris = []
    for dr in INPUT_DIR:
        sample = os.path.basename(dr)
        results_dir = os.path.join(dr, "hprc_chry/assembly_qc/nucflag/v0.3.4_hifi/")
        for file in os.listdir(results_dir):
            samples.append(sample)
            output_urls.append(os.path.join(OUTPUT_URL.format(sample=sample), file))
            output_uris.append(os.path.join(OUTPUT_URI.format(sample=sample), file))

    df = (
        pl.DataFrame({
            "sample_id": samples,
            "urls": output_urls,
            "uri": output_uris,
        })
        .with_columns(
            dtype=pl.when(pl.col("uri").str.contains("plots"))
            .then(pl.lit("plot_tarball"))
            .when(pl.col("uri").str.contains("first"))
            .then(pl.lit("bw_first"))
            .when(pl.col("uri").str.contains("second"))
            .then(pl.lit("bw_second"))
            .when(pl.col("uri").str.contains("bed"))
            .then(pl.lit("bedfile"))
        )
        .with_columns(
            dtype=pl.when(pl.col("uri").str.contains("md5"))
            .then(pl.col("dtype") + "_" + pl.lit("hash"))
            .otherwise(pl.col("dtype"))
        )
        .sort(by=["dtype", "sample_id"])
        .select("sample_id", "dtype", "urls", "uri")
    )
    df.write_csv(sys.stdout, separator="\t", include_header=True)

if __name__ == "__main__":
    raise SystemExit(main())
