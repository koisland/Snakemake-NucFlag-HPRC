import os
import sys
import glob
import polars as pl

def main():
    indir = sys.argv[1]
    outdir = sys.argv[2]

    os.makedirs(outdir, exist_ok=True)

    BASE_URI = "https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/5E9D3123-C9D3-47D4-8BB7-D82FF0DC84EA--HPRC_R2_QC_NUCFLAG/"
    SUB_URI = "/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi/"
    DTYPE_FILES = {
        "bigwig_first": "nucflag.first.bw",
        "bigwig_second": "nucflag.second.bw",
        "misassemblies": "nucflag.bed",
        "plots": "nucflag.plots.tar.gz",
    }
    exclude_samples = {
        "HG02257",
        "HG01978",
        "HG03516"
    }
    samples = [
        os.path.basename(p) for p in glob.glob(f"{indir}/*") if not os.path.basename(p) in exclude_samples
    ]
    df_metadata = pl.read_csv(
        "/project/logsdon_shared/data/HPRC/assemblies/assemblies_pre_release_v0.6.1.index.csv"
    )
    """
    {assembly_id}.nucflag.first.bw
    {assembly_id}.nucflag.second.bw
    {assembly_id}.nucflag.bed
    {assembly_id}.nucflag.plots.tar.gz
    """
    df_metadata_subset = df_metadata.filter(pl.col("sample_id").is_in(samples))

    for dtype, file_suffix in DTYPE_FILES.items():
        outfile = os.path.join(outdir, f"nucflag_hifi_{dtype}.csv")

        df_res = (
            df_metadata_subset
            .select("sample_id", "assembly_name")
            .with_columns(
                location=BASE_URI + pl.col("sample_id") + SUB_URI + pl.col("assembly_name") + f".{file_suffix}" 
            )
        )
        (
            df_res
            .write_csv(outfile)
        )


if __name__ == "__main__":
    raise SystemExit(main())
