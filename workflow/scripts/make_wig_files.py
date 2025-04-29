import os
import sys
import glob
import polars as pl

def main():
    indir = sys.argv[1]
    first_wig = sys.argv[2]
    second_wig = sys.argv[3]
    hap = sys.argv[4]

    cov_files = sorted(glob.glob(os.path.join(indir, "*.gz")))
    chroms = set(
        os.path.basename(file).split(":")[0]
        for file in cov_files
    )
    with (
        open(first_wig, mode="a") as first_fh,
        open(second_wig, mode="a") as second_fh,
    ):
        for chrom in chroms:
            chrom_hap = chrom.split("#")[1]
            if chrom_hap != hap:
                continue
            header = f"variableStep chrom={chrom}\n"
            df = (
                pl.read_csv(
                    os.path.join(indir, f"{chrom}*.gz"),
                    separator="\t",
                    has_header=True
                )
                .select("position", "first", "second")
                .with_columns(
                    pl.col("position") + 1
                )
                .sort(by="position")
            )
            first_fh.write(header)
            (
                df
                .select("position", "first")
                .write_csv(
                    first_fh, include_header=False, separator="\t"
                )
            )
            second_fh.write(header)
            (
                df
                .select("position", "second")
                .write_csv(
                    second_fh, include_header=False, separator="\t"
                )
            )

if __name__ == "__main__":
    raise SystemExit(main())
