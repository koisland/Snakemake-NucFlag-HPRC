# NucFlag for HPRC chrY
WIP. Only works on UPenn's LPC.

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

# Create first batch
python get_next_batch.py 1 > batch_1.txt

# Generate config
python generate_config.py 1 > config_batch_1.yaml

pushd ..

# Run on samples.
snakemake -np --configfile config_batch_1.yaml -j 100 --workflow-profile ~/profiles/lpc

popd

# Then create next batch
python get_next_batch.py 2 > batch_2.txt
```
