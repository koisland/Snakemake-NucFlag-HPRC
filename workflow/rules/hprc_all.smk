import os
import glob
from os.path import join
import polars as pl

VERSION = "v0.3.3"
DTYPE = "hifi"
OUTDIR = "/project/logsdon_shared/projects/HPRC/CenMAP/nucflag_hprc/results_hprc"
# Haplotypes switched. BAM is from assembly v1.0.1
# See https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables#known-issues
EXCLUDE_SAMPLES = {
    "HG02257",
    "HG01978",
    "HG03516"
}
SAMPLE_BAMS = {}
for d in glob.glob("data/aln/*"):
    if not os.path.isdir(d) or d.startswith("."):
        continue
    sm = os.path.basename(d)
    if sm in EXCLUDE_SAMPLES:
        continue
    # Sort first to always take deep concensus bams.
    bam_files = glob.glob(join(d, "*.bam"))
    bam_files.sort()

    try:    
        bam_file = next(iter(bam_files))
    except Exception:
        continue
    
    SAMPLE_BAMS[sm] = bam_file

CHROM_ALIASES = glob.glob("/project/logsdon_shared/data/HPRC/annotations/*.chromAlias.txt")
df_metadata = pl.read_csv("/project/logsdon_shared/data/HPRC/assemblies/assemblies_pre_release_v0.6.1.index.csv")
SM_NAME_KEY = dict(
    df_metadata
    .select("sample_id", "assembly_name")
    .iter_rows()
)
NAME_ASM_KEY = dict(
    df_metadata
    .select("assembly_name", "assembly_method")
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

rule convert_nucflag_to_bed9:
    input:
        mtypes="mtype.json",
        bed=rules.run_nucflag.output.misassemblies
    output:
        bed=join(
            OUTDIR,
            "{sm}_nucflag.bed",
        )
    run:
        import json
        with open(input.mtypes) as fh:
            MTYPES = json.load(fh)
        
        with (
            open(input.bed) as fh,
            open(output.bed, "wt") as ofh
        ):
            for line in fh:
                chrom, st, end, mtype = line.strip().split()
                st, end = int(st) + 1, int(end) + 1
                st, end = str(st), str(end)
                orow = [
                    chrom,
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


def get_asm_dir(wc):
    ASM_NAME = SM_NAME_KEY[wc.sm]
    METHOD = NAME_ASM_KEY[ASM_NAME]
    return f"/project/logsdon_shared/data/HPRC/assemblies/{METHOD}/{wc.sm}/"

rule make_chrom_size:
    input:
        bam_file=lambda wc: SAMPLE_BAMS[wc.sm],
        asm_dir=get_asm_dir
    output:
        first_wig=temp(join(OUTDIR, "{sm}_chrom_sizes.tsv")),
    params:
        bam_file=lambda wc, input: os.path.realpath(input.bam_file)
    shell:
        """
        cut -f 1,2 {input.asm_dir}/*.fai | \
        grep -f <(samtools view -H {params.bam_file} | grep SQ | cut -f 2 | sed 's/SN://g') > {output}
        """

rule split_cov_to_wigs:
    input:
        script="make_wig_files.py",
        cov_dir=rules.run_nucflag.output.plot_dir
    output:
        first_wig=temp(join(OUTDIR, "{sm}_first.wig")),
        second_wig=temp(join(OUTDIR, "{sm}_second.wig")),
    conda:
        "general"
    shell:
        """
        python {input.script} \
        {input.cov_dir} \
        {output.first_wig} \
        {output.second_wig}
        """

rule convert_first_wig_to_bigwig:
    input:
        wig=rules.split_cov_to_wigs.output.first_wig,
        chrom_sizes=rules.make_chrom_size.output,
    output:
        bw=join(OUTDIR, "{sm}_nucflag_first.bw"),
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
        bw=join(OUTDIR, "{sm}_nucflag_second.bw"),
        

"""
upload_folder \
    {sample_id}/hprc_r2/assembly_qc/nucflag \
        {assembly_id}.nucflag.first.bw
        {assembly_id}.nucflag.second.bw
        {assembly_id}.nucflag.bed
        {assembly_id}.nucflag.plots.tar.gz

So sample ID would be HG000099 and asssembly id would be HG00099_hap1_hprc_v1.0.1 (from the assembly index file; in the assembly id column)
"""

rule create_tarballs:
    input:
        first_bw = rules.convert_first_wig_to_bigwig.output,
        second_bw = rules.convert_second_wig_to_bigwig.output,
        bedfile = rules.convert_nucflag_to_bed9.output,
        plot_dir = rules.run_nucflag.output.plot_dir
    output:
        join(OUTDIR, "{sm}.done")
    params:
        assembly_id=lambda wc: SM_NAME_KEY[wc.sm],
        # hprc_r2/assembly_qc/nucflag/
        tmp_dir=lambda wc: join(OUTDIR, "final", wc.sm, "hprc_r2", "assembly_qc", "nucflag", f"{VERSION}_{DTYPE}")
    conda:
        "general"
    threads:
        4
    shell:
        """
        mkdir -p {params.tmp_dir}
        ln -sf {input.first_bw} {params.tmp_dir}/{params.assembly_id}.nucflag.first.bw
        ln -sf {input.second_bw} {params.tmp_dir}/{params.assembly_id}.nucflag.second.bw
        ln -sf {input.bedfile} {params.tmp_dir}/{params.assembly_id}.nucflag.bed
        cd {input.plot_dir}
        tar -czf {params.tmp_dir}/{params.assembly_id}.nucflag.plots.tar.gz *.png
        cd {params.tmp_dir}
        ls | parallel -j {threads} 'md5sum {{}} > {{}}.md5'
        touch {output}
        """

wildcard_constraints:
    sm="|".join(SAMPLE_BAMS)

rule all:
    input:
        expand(rules.run_nucflag.output, sm=SAMPLE_BAMS),
        expand(rules.convert_nucflag_to_bed9.output, sm=SAMPLE_BAMS),
        expand(rules.convert_first_wig_to_bigwig.output, sm=SAMPLE_BAMS),
        expand(rules.convert_second_wig_to_bigwig.output, sm=SAMPLE_BAMS),
        expand(rules.create_tarballs.output, sm=SAMPLE_BAMS),
    default_target: True
