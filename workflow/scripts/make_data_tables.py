import os
import sys
import glob
import polars as pl


# HAP_SWITCH_SAMPLES = {
#     "HG02257",
#     "HG01978",
#     "HG03516"
# }
def normalize_hap(asm_name: str) -> str:
    sm, hap, *_ = asm_name.split("_")

    if hap == "pat" or hap == "hap1":
        hap = "1"
    elif hap == "mat" or hap == "hap2":
        hap = "2"
    else:
        raise ValueError(hap)
    return hap


def main():
    indir = sys.argv[1]
    outdir = sys.argv[2]

    os.makedirs(outdir, exist_ok=True)

    BASE_URI = "s3://human-pangenomics/submissions/5E9D3123-C9D3-47D4-8BB7-D82FF0DC84EA--HPRC_R2_QC_NUCFLAG/"
    SUB_URI = "/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi/"
    DTYPE_FILES = {
        "bigwig_first": "nucflag.first.bw",
        "bigwig_second": "nucflag.second.bw",
        "misassemblies": "nucflag.bed",
        "plots": "nucflag.plots.tar.gz",
    }
    samples = [
        os.path.basename(p) for p in glob.glob(f"{indir}/*")
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
                haplotype=pl.col("assembly_name").map_elements(normalize_hap, return_dtype=pl.String),
                location=BASE_URI + pl.col("sample_id") + SUB_URI + pl.col("assembly_name") + f".{file_suffix}" 
            )
            .sort("haplotype", "sample_id")
            .select("sample_id", "haplotype", "assembly_name", "location")
        )
        (
            df_res
            .write_csv(outfile)
        )


if __name__ == "__main__":
    raise SystemExit(main())
