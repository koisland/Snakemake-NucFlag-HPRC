# Description

NucFlag misassembled regions 

# Methods

NucFlag is a tool to flag regions of the genome that are potentially misassembled based on read alignments.

It achieves this by:
1. Generating a read pileup of the first and second most common base frequencies per position in the genome.
2. Calculating the heterozygous base ratio.   
3. Calling peaks in both first and second signals.
4. Filtering peaks with the most abnormal read alignment coverage based on given thresholds.
5. Overlapping peaks to determine misassembly classifications.

# Misassemblies

MISJOIN (orange)
* Drop in first most common base coverage or a coverage gap.
* This region has minimal to no reads supporting it or is a scaffold.
* Can overlap region with secondary base coverage.

COLLAPSE (green)
* Collapse with no variants.

COLLAPSE_VAR (blue)
* Collapse with variants. Overlaps region with high secondary base coverage.

COLLAPSE_OTHER (red)
* Region with high het ratio.

HET (teal)
* Possible heterozygous (het) site.
* Determined by the het ratio, the coverage of the second most common base divided by the first most common base coverage plus the second most common base coverage.


# Credits

Keisuke K. Oshima <Keith.Oshima@Pennmedicine.upenn.edu>
