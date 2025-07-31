import os, sys, yaml
import glob

SM_TEMPLATE = {
    "name": "{sm}",
    "asm_dir": "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/data/asm_hgsvc/grouped/{sm}",
    "asm_rgx": ".*\\.gz$",
    "read_fofn": "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/data/hifi_hgsvc/{sm}.fofn",
    "config": "config/nucflag.toml",
}

def main():
    cfg_template = "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/config/config_template.yaml"
    with open(cfg_template) as fh:
        cfg = yaml.safe_load(fh)

    sample_cfg = []
    samples = [
        os.path.basename(path)
        for path in glob.glob("/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/data/asm_hgsvc/grouped/*")
        if os.path.isdir(path)
    ]
    for sm in samples:
        sm = sm.strip()
        sample_cfg.append({k: v.format(sm=sm) for k,v in SM_TEMPLATE.items()})

    cfg["samples"] = sample_cfg
    yaml.safe_dump(cfg, sys.stdout)

if __name__ == "__main__":
    raise SystemExit(main())
