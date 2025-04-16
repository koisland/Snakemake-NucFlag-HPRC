import ast
from matplotlib.axes import Axes
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt

import os
import json
import argparse
from os.path import join, basename
from typing import Any, Literal

CHROMS = [f"chr{chrom}" for chrom in (*range(1, 23), "X", "Y")]
CHROMS.append("Other")


def draw_nucflag_breakdown_plot(
    df_length: pl.DataFrame,
    x: str,
    y: str,
    outfile_prefix: str,
    unit: Literal["percent", "Mbp"] = "Mbp",
    palette: dict[str, str] | None = None,
    ylim: tuple[int, int] | None = None,
):
    fig, ax = plt.subplots(figsize=(21, 7))
    ax: Axes

    if unit == "percent":
        df_length = df_length.with_columns((pl.col(y) / pl.col(y).sum().over(x)) * 100)

    sns.histplot(
        df_length,
        x=x,
        weights=y,
        hue="mtype",
        multiple="stack",
        palette=palette,
        edgecolor="black",
        hue_order=sorted(df_length["mtype"].unique()),
        linewidth=1.0,
        ax=ax,
    )

    _ = plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    if unit == "Mbp":
        ax.set_ylabel(
            f"Cumulative {y} evaluated by NucFlag ({unit})"
        )
        ax.ticklabel_format(axis="y", useOffset=False, style="plain")
        yticks = ax.get_yticks()
        ax.set_yticks(yticks, [str(tk / 1_000_000) for tk in yticks])
    else:
        ax.set_ylabel(
            f"{y.capitalize()} evaluated by NucFlag (%)"
        )

    ax.set_xlabel(None)

    if ylim:
        ax.set_ylim(ylim)

    # Hide spines
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    sns.move_legend(
        ax,
        title="Types",
        alignment="left",
        # Set aspect ratio of handles so square.
        handlelength=1.0,
        handleheight=1.0,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        ncol=1,
    )
    fig.savefig(f"{outfile_prefix}.png", bbox_inches="tight")
    fig.savefig(f"{outfile_prefix}.pdf", bbox_inches="tight")
    plt.close()


def plot_summary(
    fais: list[str],
    misassemblies: list[str],
    chrom_aliases: list[str],
    outdir: str,
    mcolors: str
):
    with open(mcolors, "rt") as fh:
        mcolors: dict[str, Any] = json.load(fh)
        mcolors = {
            mtype: tuple(v / 255 for v in ast.literal_eval(f"({item_rgb})"))
            for mtype, item_rgb in mcolors.items()
        }

    os.makedirs(outdir, exist_ok=True)

    dfs_lengths = []
    for file in fais:
        sm = basename(file).split("-")[0]
        df = pl.read_csv(
            file,
            columns=[0, 1],
            new_columns=["chrom", "length"],
            separator="\t",
            has_header=False,
        ).with_columns(sm=pl.lit(sm))
        dfs_lengths.append(df)

    df_length: pl.DataFrame = pl.concat(dfs_lengths)
    df_length = df_length.with_columns(
        pl.col("chrom").str.extract(r".*?_chr[0-9XY]+_(.*)").fill_null(pl.col("chrom"))
    )

    df_chrom = (
        pl.concat(
            pl.read_csv(
                c_alias,
                has_header=True,
                separator="\t",
            )
            for c_alias in chrom_aliases
        )
        .with_columns(chrom_name=pl.col("ucsc").str.extract(r"(chr[0-9XY]+)"))
        .drop_nulls()
        .cast({"chrom_name": pl.Enum(CHROMS)})
        .rename({"# assembly": "chrom"})
        .select("chrom", "chrom_name")
        # Extract either sm_chr_ctg or ctg
        .with_columns(
            pl.col("chrom")
            .str.extract(r".*?_chr[0-9XY]+_(.*)")
            .fill_null(pl.col("chrom"))
        )
    )

    df_misassemblies_length = (
        pl.concat(
            pl.read_csv(
                msm,
                has_header=False,
                separator="\t",
                new_columns=["chrom", "st", "end", "mtype"],
            )
            for msm in misassemblies
        )
        .group_by(["chrom", "mtype"])
        .agg(length=(pl.col("end") - pl.col("st")).sum())
        # Extract either sm_chr_ctg or ctg
        .with_columns(
            pl.col("chrom")
            .str.extract(r".*?_chr[0-9XY]+_(.*)")
            .fill_null(pl.col("chrom"))
        )
        .with_columns(
            new_chrom=pl.col("chrom")
            .str.splitn("#", 2)
            .struct.rename_fields(["sm", "ctg"]),
        )
        .unnest("new_chrom")
        .sort(by=["sm"])
    )
    # sample
    df_sm_missassembled_length = (
        df_misassemblies_length
        .group_by(["sm", "mtype"])
        .agg(pl.col("length").sum())
    )
    df_sm_region_length = (
        df_length
        .group_by(["sm"])
        .agg(mtype=pl.lit("GOOD"), length=pl.col("length").sum())
        .join(
            df_misassemblies_length.group_by(["sm"]).agg(pl.col("length").sum()),
            on="sm",
            how="left",
        )
        .with_columns(pl.col("length") - pl.col("length_right").fill_null(0))
        .select("sm", "mtype", "length")
    )
    # chrom
    df_chrom_missassembled_length = (
        df_misassemblies_length
        .join(df_chrom, on="chrom", how="left")
        .group_by(["chrom_name", "mtype"])
        .agg(pl.col("length").sum())
        .with_columns(pl.col("chrom_name").fill_null(pl.lit("Other")))
    )
    df_chrom_region_length = (
        df_length
        .join(df_chrom, on="chrom", how="left")
        .group_by(["chrom_name"])
        .agg(mtype=pl.lit("GOOD"), length=pl.col("length").sum())
        .join(
            df_chrom_missassembled_length.group_by(["chrom_name"]).agg(
                pl.col("length").sum()
            ),
            on="chrom_name",
            how="left",
        )
        .with_columns(pl.col("length") - pl.col("length_right").fill_null(0))
        .select("chrom_name", "mtype", "length")
        .with_columns(pl.col("chrom_name").fill_null(pl.lit("Other")))
    )

    df_chrom_length = pl.concat(
        [df_chrom_region_length, df_chrom_missassembled_length]
    ).sort(by=["chrom_name", "mtype"])
    df_sm_length = pl.concat([df_sm_region_length, df_sm_missassembled_length]).sort(
        by=["sm", "mtype"],
    )

    df_sm_summary = (
        df_sm_length.pivot(index="sm", on="mtype", values="length")
        .with_columns(TOTAL=pl.sum_horizontal(*mcolors.keys()))
        .with_columns(
            **{
                m: ((pl.col(m) / pl.col("TOTAL")) * 100.0).fill_null(0.0)
                for m in mcolors
            }
        )
        .sort(by="sm")
    )
    df_sm_summary.write_csv(join(outdir, "sm_summary.csv"))

    df_long_sm_summary = (
        df_sm_summary.drop("TOTAL")
        .unpivot(index="sm", variable_name="mtype", value_name="percent")
        .sort(by=["sm", "percent"], descending=True)
    )

    df_chrom_summary = (
        df_chrom_length.pivot(index="chrom_name", on="mtype", values="length")
        .with_columns(TOTAL=pl.sum_horizontal(*mcolors.keys()))
        .with_columns(
            **{
                m: ((pl.col(m) / pl.col("TOTAL")) * 100.0).fill_null(0.0)
                for m in mcolors
            }
        )
        .sort(by="chrom_name")
    )
    df_chrom_summary.write_csv(join(outdir, "chrom_summary.csv"))

    df_long_chrom_summary = (
        df_chrom_summary.drop("TOTAL")
        .unpivot(index="chrom_name", variable_name="mtype", value_name="percent")
        .sort(by=["chrom_name", "percent"], descending=True)
    )

    # sample
    draw_nucflag_breakdown_plot(
        df_sm_length, "sm", "length", join(outdir, "sm_nucflag"),
        palette=mcolors
    )
    draw_nucflag_breakdown_plot(
        df_sm_length.filter(pl.col("mtype") != "GOOD"),
        "sm",
        "length",
        join(outdir, "sm_nucflag_misassemblies"),
        palette=mcolors,
        ylim=(0, 400_000_000)
    )
    draw_nucflag_breakdown_plot(
        df_long_sm_summary,
        "sm",
        "percent",
        join(outdir, "sm_nucflag_percent"),
        unit="percent",
        palette=mcolors
    )
    draw_nucflag_breakdown_plot(
        df_long_sm_summary.filter(pl.col("mtype") != "GOOD"),
        "sm",
        "percent",
        join(outdir, "sm_nucflag_misassemblies_percent"),
        unit="percent",
        palette=mcolors
    )
    df_long_sm_summary.write_csv(join(outdir, "sm_nucflag_summary.csv"))
    df_sm_length.write_csv(join(outdir, "sm_nucflag.csv"))

    # Chrom
    draw_nucflag_breakdown_plot(
        df_chrom_length, "chrom_name", "length", join(outdir, "chrom_nucflag"), palette=mcolors
    )
    draw_nucflag_breakdown_plot(
        df_chrom_length.filter(pl.col("mtype") != "GOOD"),
        "chrom_name",
        "length",
        join(outdir, "chrom_nucflag_misassemblies"),
        palette=mcolors,
        ylim=(0, 80_000_000)
    )
    draw_nucflag_breakdown_plot(
        df_long_chrom_summary,
        "chrom_name",
        "percent",
        join(outdir, "chrom_nucflag_percent"),
        unit="percent",
        palette=mcolors
    )
    draw_nucflag_breakdown_plot(
        df_long_chrom_summary.filter(pl.col("mtype") != "GOOD"),
        "chrom_name",
        "percent",
        join(outdir, "chrom_nucflag_misassemblies_percent"),
        unit="percent",
        palette=mcolors
    )
    df_long_chrom_summary.write_csv(join(outdir, "chrom_nucflag_summary.csv"))
    df_chrom_length.write_csv(join(outdir, "chrom_nucflag.csv"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fai", nargs="+")
    ap.add_argument("-m", "--misassemblies", nargs="+")
    ap.add_argument("-c", "--chrom", nargs="+")
    ap.add_argument("-o", "--out_prefix")
    ap.add_argument("--mcolors")
    args = ap.parse_args()

    """
    python summarize_misassemblies.py \
    -f /project/logsdon_shared/projects/HPRC/CenMAP/results_hifiasm_b2-3_v0.4.2/1-concat_asm/*asm-comb-dedup.fa.fai \
    -m /project/logsdon_shared/projects/HPRC/CenMAP/nucflag_hprc/results/*_misassemblies.bed \
    -c /project/logsdon_shared/data/HPRC/annotations/*Alias.txt \
    -o /project/logsdon_shared/projects/HPRC/CenMAP/nucflag_hprc/plots \
    --mcolors /project/logsdon_shared/projects/HPRC/CenMAP/nucflag_hprc/mtype.json
    """
    plot_summary(
        fais=args.fai,
        misassemblies=args.misassemblies,
        chrom_aliases=args.chrom,
        outdir=args.out_prefix,
        mcolors=args.mcolors
    )
