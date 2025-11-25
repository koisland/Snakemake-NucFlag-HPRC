
glob_fa = "/project/logsdon_shared/projects/HPRC/CenMAP_chrY/results_{grp}/2-concat_asm/{sm}-asm-renamed-reort.fa"
glob_bam = "/project/logsdon_shared/projects/HPRC/CenMAP_chrY/results_{grp}/8-cdr_finder/aln/{sm}.bam"

wcs = glob_wildcards(glob_bam)

wildcard_constraints:
    grp="|".join(wcs.grp),
    sm="|".join(wcs.sm),

rule subset_bams:
    input:
        bam=glob_bam,
        fai=f"{glob_fa}.fai"
    output:
        bam=temp("results_{grp}/{sm}_chrY_renamed.bam"),
        bai=temp("results_{grp}/{sm}_chrY_renamed.bam.bai"),
    threads:
        12
    shell:
        """
        samtools view \
        -b {input.bam} \
        --regions-file <(grep chrY {input.fai} | awk -v OFS="\\t" '{{ print $1, 0, $2 }}') \
        -o {output.bam} \
        -@ {threads}
        samtools index {output.bam} -@ {threads}
        """

rule subset_asm:
    input:
        fa=glob_fa,
        fai=f"{glob_fa}.fai"
    output:
        fa=temp("results_{grp}/{sm}_chrY_renamed.fa.gz"),
        fai=temp("results_{grp}/{sm}_chrY_renamed.fa.gz.fai"),
        gzi=temp("results_{grp}/{sm}_chrY_renamed.fa.gz.gzi"),
    threads:
        12
    shell:
        """
        seqtk subseq {input.fa} <(grep chrY {input.fai} | cut -f 1) | bgzip > {output.fa}
        samtools faidx {output.fa}
        """

rule subset_bam_only_main_Y:
    input:
        bam=rules.subset_bams.output.bam,
        fai=f"{glob_fa}.fai"
    output:
        bam="results_{grp}/{sm}_chrY.bam",
        bai="results_{grp}/{sm}_chrY.bam.bai",
        # bam=temp("results_{grp}/{sm}_chrY.bam"),
        # bai=temp("results_{grp}/{sm}_chrY.bam.bai"),
    shell:
        """
        samtools view -h {input.bam} | \
        awk -v OFS="\\t" '{{
            # Write header and rename SQ
            if ($1 ~ "^@") {{
                if ($1 ~ "@SQ") {{
                    if ($2 !~ "chrY" || $2 ~ "unassigned" || $2 ~ "random") {{ next }}
                    match($2, "([^_]*?)_chrY$", sm);
                    if (!sm[1]) {{ next }}
                    $2="SN:"sm[1]"_chrY"
                }}
                print
                next
            }}
            # No unassigned or random
            if ($3 ~ "unassigned" || $3 ~ "random") {{ next }}
            match($3, "([^_]*?)_chrY$", sm);
            # Restore original chrom name
            # This should be safe so long as no contigs are reoriented.
            if (!sm[1]) {{ next }}
            $3=sm[1]"_chrY";
            print
        }}' | \
        sed -e 's/{wildcards.sm}_chrY_{wildcards.sm}_chrY/{wildcards.sm}_chrY/g' \
            -e 's/{wildcards.sm}_rc-chrY_{wildcards.sm}_chrY/{wildcards.sm}_chrY/g' | \
        samtools view -b -o {output.bam}
        samtools index {output.bam}
        """



rule all:
    input:
        expand(rules.subset_bams.output, zip, grp=wcs.grp, sm=wcs.sm),
        expand(rules.subset_asm.output, zip, grp=wcs.grp, sm=wcs.sm),
        expand(rules.subset_bam_only_main_Y.output, zip, grp=wcs.grp, sm=wcs.sm),
    default_target:
        True
