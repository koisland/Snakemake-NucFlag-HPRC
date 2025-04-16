from os.path import join

SAMPLES = set(glob_wildcards("data/aln/{sm}.hifi_minimap2_2.28.corrected.bam").sm)
FAIS = {
    "mm2": expand(
        "/project/logsdon_shared/projects/HPRC/CenMAP/results_hifiasm_b2-3_v0.4.2/1-concat_asm/{sm}-asm-comb-dedup.fa.fai",
        sm=SAMPLES
    ),
    "pbmm2": expand(
        "/project/logsdon_shared/projects/HPRC/CenMAP/results_hifiasm_b2-3_v0.4.2/1-concat_asm/{sm}-asm-renamed-reort.fa.fai",
        sm=SAMPLES
    )
}
GLOB_CHROM_ALIAS = "/project/logsdon_shared/data/HPRC/annotations/{sm}_{hap}_hprc_r2_v1.0.1.chromAlias.txt"
WC_CHROM_ALIAS = glob_wildcards(GLOB_CHROM_ALIAS)
CHROM_ALIAS_SM, CHROM_ALIAS_HAP = [], []
for sm, hap in zip(WC_CHROM_ALIAS.sm, WC_CHROM_ALIAS.hap):
    if sm in SAMPLES:
        CHROM_ALIAS_SM.append(sm)
        CHROM_ALIAS_HAP.append(hap)
CHROM_ALIASES = expand(
    GLOB_CHROM_ALIAS,
    zip,
    sm=CHROM_ALIAS_SM,
    hap=CHROM_ALIAS_HAP,
)

rule make_rm_overlay_bed:
    input:
        bed_dir=join("data", "rm", "{sm}")
    output:
        bed=join("results", "{sm}_rm.bed")
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

# TODO: Renamed centromeric contigs in RM.
rule rename_rm_overlay_bed:
    input:
        bed_dir=join("data", "rm", "{sm}"),
        bed_key="/project/logsdon_shared/projects/HPRC/CenMAP/results_hifiasm_b2-3_v0.4.2/4-ident_cen_ctgs/bed/interm/{sm}_renamed_cens.tsv"
    output:
        bed=join("results", "{sm}_rm_renamed.bed")
    shell:
        """
        find {input.bed_dir} -name "*.bed" -exec cat {{}} \; | \
        sort | \
        join -a1 -a2 - <(sort -k1 {input.bed_key}) | \
        awk -v OFS="\\t" '{{
            if ($11 == "") {{
                new_chrom = $1
            }} else {{
                new_chrom = $11
            }}
            if ($7 == "Satellite") {{
                name=$4
            }} else {{
                name=$7
            }}
            print new_chrom, $2, $3, name, "plot"
        }}' > {output}
        """

rule run_nucflag:
    input:
        bam_file="data/aln/{sm}.hifi_minimap2_2.28.corrected.bam",
        config="nucflag.toml",
        ovl_bed=rules.make_rm_overlay_bed.output,
    output:
        plot_dir=directory(join("results", "{sm}")),
        misassemblies=join(
            "results",
            "{sm}_misassemblies.bed",
        ),
    threads: 12
    resources:
        mem="60GB",
    conda:
        "general"
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
        --overlay_regions {input.ovl_bed} &> {log}
        """

use rule run_nucflag as run_nucflag_pbmm2 with:
    input:
        bam_file="/project/logsdon_shared/projects/HPRC/CenMAP/results_hifiasm_b2-3_v0.4.2/8-nucflag/{sm}.bam",
        config="nucflag.toml",
        ovl_bed=rules.rename_rm_overlay_bed.output,
    output:
        plot_dir=directory(join("results", "{sm}_pbmm2")),
        misassemblies=join(
            "results",
            "{sm}_misassemblies_pbmm2.bed",
        ),
    log:
        os.path.join("logs", "run_nucflag_{sm}_pbmm2.log"),
    benchmark:
        os.path.join("benchmark", "run_nucflag_{sm}_pbmm2.tsv")
    

rule convert_nucflag_to_bed9:
    input:
        mtypes="mtype.json",
        bed=rules.run_nucflag_pbmm2.output.misassemblies
    output:
        bed=join(
            "results",
            "{sm}_misassemblies_pbmm2_final.bed",
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

use rule convert_nucflag_to_bed9 as convert_nucflag_to_bed9_mm2 with:
    input:
        mtypes="mtype.json",
        bed=rules.run_nucflag.output.misassemblies
    output:
        bed=join(
            "results",
            "{sm}_misassemblies_final.bed",
        )

rule intersect_beds:
    input:
        pbmm2_bed=rules.run_nucflag_pbmm2.output.misassemblies,
        mm2_bed=rules.run_nucflag.output.misassemblies,
    output:
        intersection_mm2_pbmm2=join(
            "results",
            "{sm}_misassemblies_intersection_mm2_pbmm2.bed",
        ),
    shell:
        """
        sort -k1,1 -k2,2n \
            <(
                bedtools intersect \
                -a {input.mm2_bed} \
                -b <(awk -v OFS="\\t" '{{
                    match($1, ".+_chr[XY0-9]+_(.+)", chroms);
                    chrom=((chroms[1] != "") ? chroms[1] : $1);
                    $1=chrom;
                    print
                }}' {input.pbmm2_bed}) \
                -wa -wb -loj
            ) \
            <(
                bedtools intersect \
                -a <(awk -v OFS="\\t" '{{
                    match($1, ".+_chr[XY0-9]+_(.+)", chroms);
                    chrom=((chroms[1] != "") ? chroms[1] : $1);
                    $1=chrom;
                    print
                }}' {input.pbmm2_bed}) \
                -b {input.mm2_bed} \
                -wa -wb -loj | \
                awk -v OFS="\\t" '{{print $5, $6, $7, $8, $1, $2, $3, $4}}'
            ) | \
        uniq > {output.intersection_mm2_pbmm2}
        """

rule summarize_diff:
    input:
        bed=rules.intersect_beds.output
    output:
        directory("plots_diff/{sm}")
    params:
        bed=lambda wc, input: input.bed
    conda:
        "general"
    shell:
        """
        python summarize_diff.py -i {params.bed} -o {output}
        """

use rule summarize_diff as summarize_diff_all with:
    input:
        bed=expand(rules.intersect_beds.output, sm=SAMPLES)
    output:
        directory("plots_diff/all")
    params:
        bed=lambda wc, input: f"<(cat {input.bed})"

MISASSEMBLIES = {
    "mm2": expand(
        rules.run_nucflag.output.misassemblies,
        sm=SAMPLES
    ),
    "pbmm2": expand(
        rules.run_nucflag_pbmm2.output.misassemblies,
        sm=SAMPLES
    ),
}

rule summarize_misassemblies:
    input:
        misassemblies=lambda wc: MISASSEMBLIES[wc.aln],
        fais=lambda wc: FAIS[wc.aln],
        chrom_aliases=CHROM_ALIASES,
        mcolors="/project/logsdon_shared/projects/HPRC/CenMAP/nucflag_hprc/mtype.json"
    output:
        directory("plots_summary/{aln}")
    conda:
        "general"
    shell:
        """
        python summarize_misassemblies.py \
        -f {input.fais} \
        -m {input.misassemblies} \
        -c {input.chrom_aliases} \
        -o {output} \
        --mcolors {input.mcolors}
        """

rule all:
    input:
        expand(rules.run_nucflag.output, sm=SAMPLES),
        expand(rules.run_nucflag_pbmm2.output, sm=SAMPLES),
        expand(rules.convert_nucflag_to_bed9.output, sm=SAMPLES),
        expand(rules.convert_nucflag_to_bed9_mm2.output, sm=SAMPLES),
        expand(rules.intersect_beds.output, sm=SAMPLES),
        expand(rules.summarize_diff.output, sm=SAMPLES),
        rules.summarize_diff_all.output,
        expand(rules.summarize_misassemblies.output, aln=["mm2", "pbmm2"]),
    default_target: True
