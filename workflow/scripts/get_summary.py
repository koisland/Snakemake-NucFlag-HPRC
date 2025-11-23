import os
import sys
import glob
import intervaltree as it
import polars as pl
from collections import defaultdict

WD = os.path.dirname(__file__)


def main():
    glob_calls = glob.glob(
        "/project/logsdon_shared/projects/HPRC/CenMAP_chrY/KO_working/nucflag_v1.0/final/*/hprc_chry/assembly_qc/nucflag/v1.0.0_ont/*.bed"
    )
    itrees = defaultdict(it.IntervalTree)
    for call in glob_calls:
        if "renamed" in call:
            continue
        print(f"On {call}")
        df_call = pl.read_csv(
            call,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            columns=list(range(0, 4)),
            new_columns=["#chrom", "chromStart", "chromEnd", "name"],
        )
        for row in df_call.iter_rows(named=True):
            itrees[row["#chrom"]].add(
                interval=it.Interval(row["chromStart"], row["chromEnd"], row["name"])
            )

    df_cytobands = (
        pl.read_csv(
            "/project/logsdon_shared/projects/Keith/NucFlag/exp/HPRC_chrY/cytobands.bed",
            separator="\t",
            new_columns=["#chrom", "chromStart", "chromEnd", "name", "categ"],
            has_header=False,
        )
        .filter(~pl.col("#chrom").str.contains("random"))
        .sort(by=["#chrom", "chromStart"])
        .with_columns(gid=pl.col("chromStart").rle_id().over("#chrom"))
        .with_columns(
            categ=pl.when(pl.col("gid") == 0)
            .then(pl.lit("parm"))
            .when((pl.col("gid") == 1) | (pl.col("gid") == 2))
            .then(pl.lit("acen"))
            .otherwise(pl.lit("qarm"))
        )
        .group_by(["#chrom", "categ"])
        .agg(
            pl.col("chromStart").min(),
            pl.col("chromEnd").max(),
        )
        .sort(by=["#chrom", "chromStart"])
    )

    rows = []
    for row in df_cytobands.iter_rows():
        chrom, arm, st, end = row
        length = end - st
        ovl: set[it.Interval] = itrees[chrom].overlap(st, end)
        if not ovl:
            print(f"Skipping {row}", file=sys.stderr)
            continue

        total_correct = 0
        for itv in ovl:
            if itv.data != "correct":
                continue
            # Clip to bounds of cytoband
            ovl_st = max(itv.begin, st)
            ovl_end = min(itv.end, end)
            ovl_length = ovl_end - ovl_st
            total_correct += ovl_length
        perc_correct = (total_correct / length) * 100
        rows.append([chrom, st, end, arm, perc_correct])

    df_categs = pl.DataFrame(
        rows,
        orient="row",
        schema=["#chrom", "chromStart", "chromEnd", "categ", "perc_correct"],
    )
    df_categs.write_csv(
        os.path.join(WD, "perc_correct_by_categ.tsv"),
        separator="\t",
        include_header=True,
    )

    df_categs.group_by(["categ"]).agg(pl.col("perc_correct").mean()).write_csv(
        os.path.join(WD, "all_perc_correct_by_categ.tsv"),
        separator="\t",
        include_header=True,
    )

    # df_status = pl.read_csv(
    #     "/project/logsdon_shared/projects/Keith/NucFlag/exp/HPRC_chrY/status.bed",
    #     separator="\t",
    #     comment_prefix="#",
    #     has_header=False,
    #     columns=list(range(0,5)),
    #     new_columns=["#chrom","chromStart", "chromEnd", "status", "correct"]
    # )
    # df_status = df_status.filter(~pl.col("#chrom").str.contains("random"))
    # mean_correct = df_status["correct"].mean()
    # print(f"{mean_correct=}")


if __name__ == "__main__":
    raise SystemExit(main())
