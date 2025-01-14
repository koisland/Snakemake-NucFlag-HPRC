
rule download_repeatmasker:
    output:
        directory(config["output_rm"])
    params:
        manifest=config["manifest_rm"],
    threads:
        50
    conda:
        "../env/data.yaml"
    shell:
        """
        wget -O - {params.manifest} | \
        cut -f 1,3 -d, | \
        tail -n+1 | \
        parallel -j {threads} --colsep ',' \
        'mkdir -p {output}/{{1}}/ && aws s3 --no-sign-request cp {{2}} {output}/{{1}}/'
        """

rule download_assemblies:
    output:
        directory(config["output_asm"])
    params:
        manifest=config["manifest_asm"]
    threads:
        50
    conda:
        "../env/data.yaml"
    shell:
        """
        wget -O - {params.manifest} | \
        cut -f 1,4,13 -d, | \
        tail -n+1 | \
        parallel -j {threads} --colsep ',' \
        'mkdir -p {output}/{{2}}/{{1}}/ && aws s3 --no-sign-request cp {{3}} {output}/{{2}}/{{1}}/'
        """

rule download_hifi:
    output:
        directory(config["output_hifi"])
    params:
        manifest=config["manifest_hifi"]
    threads:
        50
    conda:
        "../env/data.yaml"
    shell:
        """
        wget -O - {params.manifest} | \
        cut -f 2,3 -d, | \
        tail -n+1 | \
        parallel -j {threads} --colsep ',' \
        'mkdir -p {output}/{{1}}/ && aws s3 --no-sign-request cp {{2}} {output}/{{1}}/'
        """

checkpoint download_all:
    input:
        asm=rules.download_assemblies.output,
        hifi=rules.download_hifi.output,
        rm=rules.download_repeatmasker.output,
    output:
        touch("download_data.done")
