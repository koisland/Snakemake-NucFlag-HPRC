import os
import sys
import glob
import numpy as np
import polars as pl
import pyBigWig as pbw
from collections import defaultdict

def main():
    indir = sys.argv[1]
    first_bw = sys.argv[2]
    second_bw = sys.argv[3]

    cov_files = sorted(glob.glob(os.path.join(indir, "*.gz")))

    chrom_sizes = defaultdict(lambda: 0)
    for file in cov_files:
        chrom, coords = os.path.basename(file).replace(".bed.gz", "").split(":")
        st, end = coords.split("-")
        end = int(end)
        chrom_sizes[chrom] = max(end, chrom_sizes[chrom]) 
        
    first_bw_fh = pbw.open(first_bw, "w")
    second_bw_fh = pbw.open(second_bw, "w")

    chrom_sizes = list(kv for kv in chrom_sizes.items())
    first_bw_fh.addHeader(chrom_sizes)
    second_bw_fh.addHeader(chrom_sizes)

    for chrom, _ in chrom_sizes:
        df = (
            pl.read_csv(
                os.path.join(indir, f"{chrom}*.gz"),
                separator="\t",
                has_header=True,
                glob=True,
            )
            .select("position", "first", "second")
            .with_columns(
                pl.col("position") + 1
            )
            .sort(by="position")
        )
        starts = df["position"].to_numpy()
        ends = starts
        first = df["first"].to_numpy().astype(np.float32)
        second = df["second"].to_numpy().astype(np.float32)
        first_bw_fh.addEntries(chrom, starts, ends=ends, values=first, span=1)
        second_bw_fh.addEntries(chrom, starts, ends=ends, values=second, span=1)
    
    first_bw_fh.close()
    second_bw_fh.close()

if __name__ == "__main__":
    raise SystemExit(main())