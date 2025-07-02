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
