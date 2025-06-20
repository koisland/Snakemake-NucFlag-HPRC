import os, sys
import glob
import random

def main():
    wd = os.path.dirname(__file__)
    
    n = int(sys.argv[1])

    batches = glob.glob(os.path.join(wd, "batch_*.txt"))
    samples_done = set()
    for batch in batches:
        with open(batch, "rt") as fh:
            for line in fh:
                line = line.strip()
                samples_done.add(line)

    all_samples = set()
    with open(os.path.join(wd, "samples.txt")) as fh:
        for line in fh:
            line = line.strip()
            all_samples.add(line)

    remaining_samples = list(all_samples.difference(samples_done)  )

    next_batch = random.sample(remaining_samples, n)

    for sm in next_batch:
        print(sm)

if __name__ == "__main__":
    raise SystemExit(main())
