import os, sys, glob


def main():
    fa_chrY_dir = "/project/logsdon_shared/data/HPRC/assemblies/verkko_chrY/chrY"
    fa_non_chrY_dir = "/project/logsdon_shared/data/HPRC/assemblies/verkko_chrY/non_chrY"
    samples = "/project/logsdon_shared/data/HPRC/assemblies/verkko_chrY/samples.txt"
    outdir = "/project/logsdon_shared/data/HPRC/assemblies/verkko_chrY/grouped"
    os.makedirs(outdir, exist_ok=True)

    with open(samples) as fh:
        samples = [line.strip() for line in fh] 
    
    for sm in samples:
        new_dir = os.path.join(outdir, sm)
        os.makedirs(new_dir, exist_ok=True)

        fa_chrY = os.path.join(fa_chrY_dir, f"{sm}_chrY.fa.gz")
        dest_fa_chrY = os.path.join(new_dir, f"{sm}_chrY.fa.gz")
        
        fa_non_chrY = os.path.join(fa_non_chrY_dir, f"{sm}_nonY.fa.gz")
        dest_fa_non_chrY = os.path.join(new_dir, f"{sm}_nonY.fa.gz")

        os.symlink(fa_chrY, dest_fa_chrY)
        os.symlink(fa_non_chrY, dest_fa_non_chrY)


if __name__ == "__main__":
    raise SystemExit(main())
