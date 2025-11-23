# NucFlag for HPRC chrY
Configured for UPenn's LPC. Replace paths in related scripts and workflow as needed.

## Output
See `results.tsv`

## Setup & Usage
Run workflow in batches.
```bash
# chrY
aws s3 sync --no-sign-request s3://human-pangenomics/T2T/scratch/chrY/v2/chrY_assemblies/ /project/logsdon_shared/data/HPRC/assemblies/verkko_chrY/chrY
aw s3 sync --no-sign-request s3://human-pangenomics/T2T/scratch/chrY/v2/nonY_assemblies/ /project/logsdon_shared/data/HPRC/assemblies/verkko_chrY/non_chrY

# Generate samples list.
find /project/logsdon_shared/data/HPRC/assemblies/verkko_chrY/chrY -name "*.fa.gz" -exec basename {} _chrY.fa.gz \; > config/samples.txt

# Then symlink files.
python config/symlink_to_grouped.py
ln -s /project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/data data/hifi
ln -s /project/logsdon_shared/data/HPRC/assemblies/verkko_chrY/grouped data/asm

# Generate configfiles
pushd config

# Create first batch of 25.
python get_next_batch.py 25 > batch_1.txt

# Generate config for batch 1
python generate_config.py 1 > config_batch_1.yaml

pushd ..

# Run on samples.
snakemake -np --configfile config_batch_1.yaml -j 100 --workflow-profile ~/profiles/lpc

popd

# Then create next batch and repeat.
python get_next_batch.py 25 > batch_2.txt

# ...
```

Upload to HPRC AWS buckets. Requires AWS credentials and submission id.
```bash
pip install git+https://github.com/DataBiosphere/ssds
bash upload_to_hprc.sh
```

Make `results.tsv`.
```bash
pip install polars
python make_data_manifest.py
```

# NucFlag ONT for HPRC chrY
Using v1.0 of NucFlag

## Run CenMAP v1.1.1
Update config with needed data. We need to run ONT alignment and use results for CDR-Finder and NucFlag.
```bash
conda create --name cenmap bioconda::cenmap
conda activate cenmap
# Then generate config.
cenmap --generate-config > hprc_chry.yaml
cenmap --generate-config > ceph_chry.yaml
cenmap --generate-config > hgsvc_chry.yaml
# Update with paths to data.
```
Run until completion.
```bash
cenmap -c $YAML -j 50 --workflow-profile ~/profiles/lpc_all
```

## Then run separate workflow for NucFlag v1.0.
```bash
# Installs via conda
snakemake -np -s workflow/rules/run_nf_v1.0.smk -j 50 --workflow-profile ~/profiles/lpc
```

## Generate ideogram
First generate cytobands from CenMAP results.
```bash
which bedtools
bash workflow/scripts/hprc_generate_cytobands.sh
```

Then generate ideogram.
```bash
which nucflag
bash workflow/scripts/generate_ideogram.sh
```

## Generate breakdown
```bash
which nucflag
bash workflow/scripts/generate_breakdown.sh
```

## Format files for upload
Format files for upload
* Fixes off-by-one in previous version.
* Converts CenMAP assembly coordinates into original coordinates.
* Computes md5 hash.
```bash
python workflow/scripts/format_for_hprc.py
```

Also had to manually fix (remove) one issue with null coordinate occuring in one sample. Recalculated md5 hash.
* `NA12886.nucflag.bed`
* TODO: Debug in NucFlag code.
```
chr14_pat-0000307	100000000	100000000	dinucleotide	0	.	100000000	100000000	0,49,83
```

## Generate summary
Uses cytobands and calls to get breakdown by region.
```bash
python /project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/workflow/scripts/get_summary.py
```
