import os, sys, yaml
import glob

SM_TEMPLATE = {
    "name": "{sm}",
    "asm_dir": "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/data/grouped_asm_ceph/{sm}",
    "asm_rgx": ".*\\.gz$",
    "read_fofn": "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/data/hifi_ceph/{sm}.fofn",
    "config": "config/nucflag.toml",
}

def main():
    wd = os.path.dirname(__file__)

    sms = glob.glob(
        "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/data/grouped_asm_ceph/*"
    )

    cfg_template = "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/config/config_template.yaml"
    with open(cfg_template) as fh:
        cfg = yaml.safe_load(fh)

    sample_cfg = []
    for sm in sms:
        sm = os.path.basename(sm)
        sample_cfg.append({k: v.format(sm=sm) for k,v in SM_TEMPLATE.items()})
    
    cfg["samples"] = sample_cfg
    yaml.safe_dump(cfg, sys.stdout)

if __name__ == "__main__":
    raise SystemExit(main())
