
glob_fa = "/project/logsdon_shared/projects/HPRC/CenMAP_chrY/results_{grp}/2-concat_asm/{sm}-asm-renamed-reort.fa"
glob_bam = "/project/logsdon_shared/projects/HPRC/CenMAP_chrY/results_{grp}/8-cdr_finder/aln/{sm}.bam"

wcs = glob_wildcards(glob_bam)

wildcard_constraints:
    grp="|".join(wcs.grp),
    sm="|".join(wcs.sm),

rule run_nucflag:
    input:
        bam=glob_bam,
        fa=glob_fa
    output:
        bed="results_{grp}/{sm}.bed",
    conda:
        "../env/nucflag_v1.0.yaml"
    params:
        preset="ont_r9"
    threads:
        12
    shell:
        """
        nucflag call -i {input.bam} -f {input.fa} -p {threads} -t {threads} -x {params.preset} -o {output.bed}
        """

rule all:
    input:
        expand(rules.run_nucflag.output, zip, grp=wcs.grp, sm=wcs.sm)
    default_target:
        True
