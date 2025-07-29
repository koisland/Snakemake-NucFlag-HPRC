import os
import glob
from os.path import join
import polars as pl

VERSION = "v0.3.3"
DTYPE = "hifi"
OUTDIR = "/project/logsdon_shared/projects/HPRC/CenMAP/nucflag_hprc/results_hprc"
TEMP_OUTDIR = join("temp", "{sm}")
FINAL_SM_OUTDIR = join(OUTDIR, "final", "{sm}", "hprc_r2", "assembly_qc", "nucflag", f"{VERSION}_{DTYPE}")
# Haplotypes switched. BAM is from assembly v1.0.1
# See https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables#known-issues
INCLUDE_SAMPLES = {
    "HG02257",
    "HG01978",
    "HG03516"
}
df_metadata = pl.read_csv("/project/logsdon_shared/data/HPRC/assemblies/assemblies_pre_release_v0.6.1.index.csv")

NAME_ASM_KEY = dict(
    df_metadata
    .select("assembly_name", "assembly_method")
    .iter_rows()
)
SAMPLE_BAMS = {}
for d in glob.glob("data/aln/*"):
    if not os.path.isdir(d) or d.startswith("."):
        continue
    sm = os.path.basename(d)
    if sm not in INCLUDE_SAMPLES:
        continue
    # Sort first to always take deep concensus bams.
    bam_files = glob.glob(join(d, "*.bam"))
    bam_files.sort()
    try:
        bam_file = next(iter(bam_files))
    except Exception:
        continue
    SAMPLE_BAMS[sm] = bam_file

SAMPLES, ASM_NAMES = zip(*
    df_metadata
    .filter(pl.col("sample_id").is_in(SAMPLE_BAMS.keys()))
    .select("sample_id", "assembly_name")
    .iter_rows()
)

rule make_rm_overlay_bed:
    input:
        bed_dir=join("data", "rm", "{sm}")
    output:
        bed=join(OUTDIR, "{sm}_rm.bed")
    shell:
        """
        find {input.bed_dir} -name "*.bed" -exec cat {{}} \; | \
        awk -v OFS="\\t" '{{
            if ($7 == "Satellite") {{
                name=$4
            }} else {{
                name=$7
            }}
            print $1, $2, $3, name, "plot"
        }}' > {output.bed}
        """

def check_rm(wc):
    if os.path.exists(join("data", "rm", f"{wc.sm}")):
        return expand(rules.make_rm_overlay_bed.output, sm=wc.sm)
    else:
        return []

rule run_nucflag:
    input:
        bam_file=lambda wc: SAMPLE_BAMS[wc.sm],
        config="nucflag.toml",
        ovl_bed=check_rm,
    output:
        plot_dir=directory(join(OUTDIR, "{sm}")),
        misassemblies=join(
            OUTDIR,
            "{sm}_misassemblies.bed",
        ),
    threads: 12
    resources:
        mem="60GB",
    conda:
        "general"
    params:
        rm_bed=lambda wc, input: f"--overlay_regions {input.ovl_bed}" if input.ovl_bed else "" 
    log:
        os.path.join("logs", "run_nucflag_{sm}.log"),
    benchmark:
        os.path.join("benchmark", "run_nucflag_{sm}.tsv")
    shell:
        """
        nucflag \
        -i {input.bam_file} \
        -d {output.plot_dir} \
        -o {output.misassemblies} \
        -t {threads} \
        -p {threads} \
        -c {input.config} \
        --output_cov_dir {output.plot_dir} \
        {params.rm_bed} &> {log}
        """

def normalize_hap(wc, rev_hap: bool):
    required_hap = wc.asm.split("_")[1]
    if required_hap == "mat" or required_hap == "hap2":
        required_hap = "1" if rev_hap else "2"
    elif required_hap == "pat" or required_hap == "hap1":
        required_hap = "2" if rev_hap else "1"
    else:
        raise ValueError(wc)
    return required_hap


def get_fai(wc):
    METHOD = NAME_ASM_KEY[wc.asm]
    return glob.glob(f"/project/logsdon_shared/data/HPRC/assemblies/{METHOD}/{wc.sm}/{wc.asm}*.fai")

rule make_chrom_size:
    input:
        # BAM has old name
        bam_file=lambda wc: SAMPLE_BAMS[wc.sm],
        # FAI has new name.
        fai=get_fai
    output:
        temp(join(TEMP_OUTDIR, "{asm}_chrom_sizes.tsv")),
    params:
        bam_file=lambda wc, input: os.path.realpath(input.bam_file)
    shell:
        """
        cut -f 1,2 {input.fai} | \
        grep -f <(samtools view -H {params.bam_file} | grep SQ | cut -f 2 | sed 's/SN://g' | awk '{{ if ($1 ~ "#1#") {{ sub("#1#", "#2#", $1) }} else {{ sub("#2#", "#1#", $1) }}; print}}') | \
        awk -v OFS="\\t" '{{
            new_name=$1;
            old_name=$1;
            if ($1 ~ "#1#") {{ sub("#1#", "#2#", new_name) }} else {{ sub("#2#", "#1#", new_name) }};
            print old_name, $2, new_name 
        }}' > {output}
        """

rule convert_nucflag_to_bed9:
    input:
        mtypes="mtype.json",
        bed=rules.run_nucflag.output.misassemblies,
        chrom_sizes=rules.make_chrom_size.output[0],
    output:
        bed=join(
            FINAL_SM_OUTDIR,
            "{asm}.nucflag.bed",
        )
    run:
        import json
        with open(input.mtypes) as fh:
            MTYPES = json.load(fh)

        rename_key = {}
        with open(input.chrom_sizes, "rt") as fh:
            for line in fh:
                old_name, _, new_name = line.strip().split("\t")
                # Need old name wig file 
                rename_key[new_name] = old_name

        with (
            open(input.bed) as fh,
            open(output.bed, "wt") as ofh
        ):
            for line in fh:
                chrom, st, end, mtype = line.strip().split()
                if chrom not in rename_key:
                    continue

                new_chrom = rename_key[chrom]

                st, end = int(st) + 1, int(end) + 1
                st, end = str(st), str(end)
                orow = [
                    new_chrom,
                    st,
                    end,
                    mtype,
                    "0",
                    ".",
                    st,
                    end,
                    MTYPES[mtype]
                ]
                ofh.write("\t".join(orow) + "\n")

rule split_cov_to_wigs:
    input:
        script="make_wig_files.py",
        cov_dir=rules.run_nucflag.output.plot_dir,
        chrom_sizes=rules.make_chrom_size.output
    output:
        first_wig=temp(join(TEMP_OUTDIR, "{asm}.nucflag.first.wig")),
        second_wig=temp(join(TEMP_OUTDIR, "{asm}.nucflag.second.wig")),
    conda:
        "general"
    shell:
        """
        python {input.script} \
        {input.cov_dir} \
        {output.first_wig} \
        {output.second_wig} \
        {input.chrom_sizes}
        """

rule convert_first_wig_to_bigwig:
    input:
        wig=rules.split_cov_to_wigs.output.first_wig,
        chrom_sizes=rules.make_chrom_size.output,
    output:
        bw=join(FINAL_SM_OUTDIR, "{asm}.nucflag.first.bw"),
    conda:
        "general"
    resources:
        mem="150GB"
    shell:
        """
        wigToBigWig {input.wig} {input.chrom_sizes} {output.bw}
        """

use rule convert_first_wig_to_bigwig as convert_second_wig_to_bigwig with:
    input:
        wig=rules.split_cov_to_wigs.output.second_wig,
        chrom_sizes=rules.make_chrom_size.output,
    output:
        bw=join(FINAL_SM_OUTDIR, "{asm}.nucflag.second.bw"),
        

"""
upload_folder \
    {sample_id}/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi \
        {assembly_id}.nucflag.first.bw
        {assembly_id}.nucflag.second.bw
        {assembly_id}.nucflag.bed
        {assembly_id}.nucflag.plots.tar.gz

So sample ID would be HG000099 and asssembly id would be HG00099_hap1_hprc_v1.0.1 (from the assembly index file; in the assembly id column)
"""

rule create_plot_tarballs:
    input:
        plot_dir = rules.run_nucflag.output.plot_dir
    output:
        join(FINAL_SM_OUTDIR, "{asm}.nucflag.plots.tar.gz")
    params:
        hap=lambda wc: normalize_hap(wc, True)
    conda:
        "general"
    shell:
        """
        cd {input.plot_dir}
        tar -czf {output} {wildcards.sm}#{params.hap}*.png
        """

rule generate_md5_hash:
    input:
        first_bw = rules.convert_first_wig_to_bigwig.output,
        second_bw = rules.convert_second_wig_to_bigwig.output,
        bedfile = rules.convert_nucflag_to_bed9.output,
        plots = rules.create_plot_tarballs.output
    output:
        rules.convert_first_wig_to_bigwig.output[0] + ".md5",
        rules.convert_second_wig_to_bigwig.output[0] + ".md5",
        rules.convert_nucflag_to_bed9.output[0] + ".md5",
        rules.create_plot_tarballs.output[0] + ".md5",
    params:
        indir = lambda wc, input: os.path.dirname(input[0])
    threads:
        4
    shell:
        """
        cd {params.indir}
        find . -mindepth 1 -not -name "*.md5" | parallel -j {threads} 'md5sum {{}} > {{}}.md5'
        """

wildcard_constraints:
    sm="|".join(SAMPLES),
    asm="|".join(ASM_NAMES),

rule all:
    input:
        expand(rules.run_nucflag.output, sm=SAMPLES),
        expand(rules.convert_nucflag_to_bed9.output, zip, sm=SAMPLES, asm=ASM_NAMES),
        expand(rules.convert_first_wig_to_bigwig.output, zip, sm=SAMPLES, asm=ASM_NAMES),
        expand(rules.convert_second_wig_to_bigwig.output, zip, sm=SAMPLES, asm=ASM_NAMES),
        expand(rules.create_plot_tarballs.output, zip, sm=SAMPLES, asm=ASM_NAMES),
        expand(rules.generate_md5_hash.output, zip, sm=SAMPLES, asm=ASM_NAMES),
    default_target: True

"""
        python make_wig_files.py         /project/logsdon_shared/projects/HPRC/CenMAP/nucflag_hprc/results_hprc/HG03516         /project/logsdon_shared/projects/HPRC/CenMAP/nucflag_hprc/results_hprc/final/HG03516/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi/HG03516_pat_hprc_r2_v1.1.0.nucflag.first.wig         /project/logsdon_shared/project
s/HPRC/CenMAP/nucflag_hprc/results_hprc/final/HG03516/hprc_r2/assembly_qc/nucflag/v0.3.3_hifi/HG03516_pat_hprc_r2_v1.1.0.nucflag.second.wig         2         true  
"""