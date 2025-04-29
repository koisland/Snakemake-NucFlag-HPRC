# Snakemake-NucFlag
[![CI](https://github.com/logsdon-lab/Snakemake-NucFlag/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/Snakemake-NucFlag/actions/workflows/main.yml)

A workflow to map PacBio HiFi reads back to an assembly and check for misassemblies at specific regions via [`NucFlag`](https://github.com/logsdon-lab/NucFlag).


### Getting Started
```bash
git clone git@github.com:logsdon-lab/Snakemake-NucFlag.git
cd Snakemake-NucFlag
```

### Configuration
Configuration is handled via `nucflag.toml`.

To read more, refer to the [`NucFlag` documentation](https://github.com/logsdon-lab/NucFlag/wiki/2.-Configuration).

#### Input
Files can be passed multiple ways in the `samples` section of `config.yaml`:

##### Assemblies
By path.
```yaml
samples:
-   name: "1"
    asm_fa: "1.fa"
```

By `fofn`.
```yaml
samples:
-   name: "1"
    asm_fofn: "1.fofn"
```

By directory and file extension.
```yaml
samples:
-   name: "1"
    asm_dir: "1/"
    asm_ext: "fa.gz"
```

##### Reads
By `fofn`.
```yaml
samples:
-   name: "1"
    read_fofn: "1.fofn"
```

By directory and file extension.
```yaml
samples:
-   name: "1"
    read_dir: "1/"
    read_ext: "bam"
```

#### Configuration

##### General
General configuration can be filled in `config.yaml`:
```yaml
# Output directory
output_dir: "results/nucflag"
# Output 1st and 2nd base coverage in {output_dir}/{sm}_coverage.
output_coverage: false
# Log directory
logs_dir: "logs/nucflag"
# Benchmarks directory
benchmarks_dir: "benchmarks/nucflag"
# Job resources. Memory in GB.
threads_aln: 8
mem_aln: 30G
processes_nucflag: 12
mem_nucflag: 50G
# samtools view filter flag.
samtools_view_flag: 2308
```

##### By Sample
The following optional are optional per sample:
```yaml
samples: [
    {
        # Regions to check. If omitted, checks entire assembly.
        region_bed: "",
        # nucflag configuration
        config: "",
        # Regions to ignore.
        ignore_bed: "",
        # Regions to overlap.
        overlay_beds: []
    }
]
```

### Output
|Path|Description|
|-|-|
|`./{output_dir}/{sample}/{contig}.png`|Per-base coverage graph plot with heterozygous sites of read coverage and potential misassemblies highlighted.|
|`./{output_dir}/{sample}_misassemblies.bed`|Bed file with heterozygous sites of read coverage and potential misassemblies with their coordinates per contig.|
|`./{output_dir}/{sample}_status.bed`|Bed file with each centromeric contig, coordinates, and status. Either `good` or `misassembled`.|
|`./{output_dir}/{sample}_coverage/{contig}.tsv`|(Optional) TSV file with base position, 1st and 2nd base coverage, and status. Either `good` or `misassembled`.|


### Usage
```bash
snakemake -np -c 1 --configfile config/config.yaml
```

### Module
To incorporate this into a workflow.

```python
SAMPLE_NAMES = ["sample_1"]
NUCFLAG_CFG = {
    "samples": [
        {
            "name": sm,
            "asm_fa": f"{sm}.fa",
            "read_dir": f"reads/{sm}/",
            "read_ext": "bam",
            "region_bed": f"regions/{sm}_region.bed",
            "overlay_beds": []
        }
        for sm in SAMPLE_NAMES
    ],
    # Other nucflag parameters
    **config["nucflag"]
}


module NucFlag:
    snakefile:
        github(
            "logsdon-lab/Snakemake-NucFlag",
            path="workflow/Snakefile",
            branch="main"
        )
    config: NUCFLAG_CFG

use rule * from NucFlag as *

rule all:
    input:
        expand(rules.nucflag.input, sm=SAMPLE_NAMES),
```
