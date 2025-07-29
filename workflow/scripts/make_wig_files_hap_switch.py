import os
import sys
import glob
import polars as pl

def main():
    indir = sys.argv[1]
    first_wig = sys.argv[2]
    second_wig = sys.argv[3]
    chrom_key = sys.argv[4]

    cov_files = sorted(glob.glob(os.path.join(indir, "*.gz")))
    chroms = set(
        os.path.basename(file).split(":")[0]
        for file in cov_files
    )

    rename_key = {}
    with open(chrom_key, "rt") as fh:
        for line in fh:
            old_name, _, new_name = line.strip().split("\t")
            # Need old name wig file 
            rename_key[new_name] = old_name

    with (
        open(first_wig, mode="a") as first_fh,
        open(second_wig, mode="a") as second_fh,
    ):
        for chrom in chroms:
            if chrom not in rename_key:
                continue

            new_chrom = rename_key[chrom]
            header = f"variableStep chrom={new_chrom}\n"
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
