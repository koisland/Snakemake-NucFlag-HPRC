import os
import argparse
from matplotlib.axes import Axes
import polars as pl
import matplotlib.pyplot as plt
import matplotlib_venn

MTYPES = (
    "HET",
    "MISJOIN",
    "COLLAPSE",
    "COLLAPSE_VAR",
    "COLLAPSE_OTHER",
)
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--infile", type=argparse.FileType("rb"))
    ap.add_argument("-o", "--output_prefix", type=str)
    args = ap.parse_args()

    os.makedirs(args.output_prefix, exist_ok=True)

    df_misassemblies_a_b = pl.read_csv(
        args.infile,
        has_header=False,
        new_columns=[
            "chrom",
            "st",
            "end",
            "mtype",
            "other_chrom",
            "other_st",
            "other_end",
            "other_mtype",
        ],
        separator="\t",
    ).with_columns(
        length=pl.col("end") - pl.col("st"),
        other_length=pl.col("other_end") - pl.col("other_st"),
    )

    df_a = df_misassemblies_a_b.filter(pl.col("other_chrom").eq("."))
    df_b = df_misassemblies_a_b.filter(pl.col("chrom").eq("."))
    # Aggregate rows with same interval on either side.
    """
    HG02293#2#CM086241.1	134008416	134015743	MISJOIN	HG02293#2#CM086241.1	134011243	134011244	MISJOIN
    HG02293#2#CM086241.1	134008416	134015743	MISJOIN	HG02293#2#CM086241.1	134011294	134011295	MISJOIN
    HG02293#2#CM086241.1	134008416	134015743	MISJOIN	HG02293#2#CM086241.1	134011309	134011310	MISJOIN
    """
    df_shared = (
        df_misassemblies_a_b.filter(
            pl.col("other_chrom").ne(".") & pl.col("chrom").ne(".")
        )
        .with_row_index()
        .with_columns(
            id=pl.col("chrom")+pl.col("st").cast(pl.String)+pl.col("end").cast(pl.String),
            other_id=pl.col("other_chrom")+pl.col("other_st").cast(pl.String)+pl.col("other_end").cast(pl.String),
        )
        # Two new cols: id - same as prev and fwd for a, other_id - same as id but for b
        .with_columns(
            id=pl.col("id").eq(pl.col("id").shift(-1)) | pl.col("id").eq(pl.col("id").shift(1)),
            other_id=pl.col("other_id").eq(pl.col("other_id").shift(-1)) | pl.col("other_id").eq(pl.col("other_id").shift(1)),
        )
        .with_columns(
            group=pl.when(
                pl.col("id") | pl.col("other_id")
            ).
            then(-1)
            .otherwise(pl.col("index"))
            .rle_id()
        )
    )
    df_shared = df_shared.group_by(["group"]).agg(
        pl.col("chrom").first(),
        pl.col("st").min(),
        pl.col("end").max(),
        # Is group the duplicated entry.
        pl.when(pl.col("id").first()).then(
            pl.col("length").first()
        ).otherwise(
            pl.col("length").sum()
        ),
        pl.col("mtype").mode().first(),
        pl.col("other_chrom").first(),
        pl.col("other_st").min(),
        pl.col("other_end").max(),
        # Is group the duplicated entry.
        pl.when(pl.col("other_id").first()).then(
            pl.col("other_length").first()
        ).otherwise(
            pl.col("other_length").sum()
        ),
        pl.col("other_mtype").mode().first(),
    )

    matplotlib_venn.venn2(
        {"10": df_a.shape[0], "01": df_b.shape[0], "11": df_shared.shape[0]},
        set_labels=["mm2", "pbmm2", "shared"],
        subset_label_formatter=lambda v: f"{v:,}",
    )
    
    plt.savefig(os.path.join(args.output_prefix, "events.png"), bbox_inches="tight")
    plt.close()
    matplotlib_venn.venn2(
        {
            "10": df_a["length"].sum(),
            "01": df_b["other_length"].sum(),
            "11": int(
                ((df_shared["length"] + df_shared["other_length"]) / 2).round().sum()
            ),
        },
        subset_label_formatter=lambda v: f"{v:,}",
        set_labels=["mm2", "pbmm2", "shared"],
    )
    plt.savefig(os.path.join(args.output_prefix, "length.png"), bbox_inches="tight")
    plt.close()

    fig, axes = plt.subplots(nrows=1, ncols=len(MTYPES), figsize=(16, 8))

    for i, mtype in enumerate(MTYPES):
        df_a_mtype = df_a.filter(pl.col("mtype") == mtype)
        df_b_mtype = df_b.filter(pl.col("other_mtype") == mtype)
        expr_all_mtype = pl.col("mtype").eq(mtype) & pl.col("other_mtype").eq(mtype)
        df_shared_mtype = df_shared.filter(expr_all_mtype)
        df_not_shared_mtype = df_shared.filter(~expr_all_mtype)
        df_a_not_shared_mtype = df_not_shared_mtype.filter(pl.col("mtype") == mtype)
        df_b_not_shared_mtype = df_not_shared_mtype.filter(pl.col("other_mtype") == mtype)
        ax: Axes = axes[i]
        venn = matplotlib_venn.venn2(
            {
                "10": df_a_mtype["length"].sum() + df_a_not_shared_mtype["length"].sum(),
                "01": df_b_mtype["other_length"].sum() + df_b_not_shared_mtype["length"].sum(),
                "11": int(
                    ((df_shared_mtype["length"] + df_shared_mtype["other_length"]) / 2)
                    .round()
                    .sum()
                ),
            },
            subset_label_formatter=lambda v: f"{v:,}",
            set_labels=["mm2", "pbmm2", "shared"],
            ax=ax
        )
        for text in venn.subset_labels:
            text.set_fontsize("x-small")
        ax.set_title(mtype)

    fig.savefig(os.path.join(args.output_prefix, f"length_by_mtype.png"), bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    raise SystemExit(main())
