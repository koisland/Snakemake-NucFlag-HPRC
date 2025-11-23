import os
import glob
import subprocess
import polars as pl

wd = os.path.dirname(__file__)
# s3://human-pangenomics/submissions/283198D5-BB8F-49CB-B1F9-CE1500812D8E--HPRC_R2_CHRY_QC_NUCFLAG/HG00096/hprc_chry/assembly_qc/nucflag/v0.3.4_hifi/
OUTPUT_DIR = "final/{sm}/hprc_chry/assembly_qc/nucflag/{version}_{dtype}/"
OUTPUT_FNAME = "{sm}.nucflag.bed"
OUTPUT_FNAME_RENAMED = "{sm}.nucflag.renamed.bed"

def main():
    linked_files = []

    # Get fasta files
    dfs_fai = {}
    for file in glob.glob(
        "/project/logsdon_shared/projects/HPRC/CenMAP_chrY/results_*/2-concat_asm/*-asm-renamed-reort.fa.fai"
    ):
        fname = os.path.basename(file)
        sm, *_ = fname.split("-")
        dfs_fai[sm] = pl.read_csv(file, separator="\t", new_columns=["#chrom", "length"], columns=[0, 1], has_header=False)

    for file in glob.glob(os.path.join(wd, "results_*", "*.bed")):
        sm, _ = os.path.splitext(os.path.basename(file))
        output_dir = OUTPUT_DIR.format(sm=sm, version="v1.0.0", dtype="ont")
        outfile = os.path.join(output_dir, OUTPUT_FNAME.format(sm=sm))
        outfile_renamed = os.path.join(output_dir, OUTPUT_FNAME_RENAMED.format(sm=sm))
        os.makedirs(output_dir, exist_ok=True)

        # Convert to original coordinate system
        df_calls = pl.read_csv(file, separator="\t", has_header=True)
        df_fai = dfs_fai[sm]

        # Fix final interval being 1 off.
        # This was fixed in https://github.com/logsdon-lab/rs-nucflag/commit/16319e885910db0146c0d9ce6605b62eb4dd4830
        df_calls = df_calls.with_columns(
            chromStart=pl.when((pl.col("chromStart") != pl.col("chromEnd").shift(1)).over("#chrom"))
            .then(pl.col("chromEnd").shift(1))
            .otherwise(pl.col("chromStart"))
        )
        df_calls_og_coords = (
            df_calls
            .join(df_fai, on="#chrom", how="left")
            .with_columns(
                is_rc=pl.col("#chrom").str.contains("rc-chr"),
                **{"#chrom": pl.col("#chrom").str.extract(r"^.*?_(rc-)*chr.*?_(.+)$", group_index=2).fill_null(pl.col("#chrom"))},
            )
            .with_columns(
                chromStart=pl.when(pl.col("is_rc"))
                .then(
                    pl.col("length")-pl.col("chromEnd")
                )
                .otherwise(pl.col("chromStart")),
                chromEnd=pl.when(pl.col("is_rc"))
                .then(
                    pl.col("length")-pl.col("chromStart")
                )
                .otherwise(pl.col("chromEnd")),
            )
            .with_columns(
                thickStart=pl.col("chromStart"),
                thickEnd=pl.col("chromEnd"),
            )
            .drop("is_rc", "length")
            .sort(by=["#chrom", "chromStart"])
        )
        assert df_calls.shape == df_calls_og_coords.shape

        df_calls.write_csv(outfile_renamed, separator="\t", include_header=True)
        df_calls_og_coords.write_csv(outfile, separator="\t", include_header=True)
        linked_files.append(outfile)
        linked_files.append(outfile_renamed)
    
    # Then compute md5sum
    for link in linked_files:
        link_dir = os.path.dirname(link)
        fname = os.path.basename(link)
        md5sum = subprocess.run(["md5sum", fname], cwd=link_dir, check=True, capture_output=True, text=True)
        with open(f"{link}.md5", "wt") as fh:
            fh.write(md5sum.stdout)

if __name__ == "__main__":
    raise SystemExit(main())
