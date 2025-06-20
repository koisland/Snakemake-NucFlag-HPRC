import os, sys, yaml
import glob

SM_TEMPLATE = {
    "name": "{sm}",
    "asm_dir": "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/data/asm/{sm}",
    "asm_rgx": ".*\\.gz$",
    "read_fofn": "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/data/hifi/{sm}.fofn",
    "config": "config/nucflag.toml",
}

def main():
    wd = os.path.dirname(__file__)
    batch_n = sys.argv[1]

    sm_batches = glob.glob(os.path.join(wd, "batches", f"batch_{batch_n}.txt"))

    cfg_template = "/project/logsdon_shared/projects/HPRC/Snakemake-NucFlag-HPRC-chrY/config/config.yaml"
    with open(cfg_template) as fh:
        cfg = yaml.safe_load(fh)

    sample_cfg = []
    for batch in sm_batches:
        with open(batch) as fh:
            for sm in fh:
                sm = sm.strip()
                sample_cfg.append({k: v.format(sm=sm) for k,v in SM_TEMPLATE.items()})
    
    cfg["samples"] = sample_cfg
    yaml.safe_dump(cfg, sys.stdout)

if __name__ == "__main__":
    raise SystemExit(main())
