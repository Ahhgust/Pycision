# pycision (V 0.001)

## What it does

pycision takes a bed file (PrecisionID_mtDNA_WG_targets.bed), 1 or more bam files and for each bam it: <br>
- infers which reads belong to which target interval in the bed <br>
  - Reads with mapping quality < Q can be optionally excluded (option: -q Q)
- trims the reads to the intervals in the bed file <br>	  
  - trimming by default only retains reads that span an entire interval <br>
  - can optionally keep reads that span at least 50% of an interval (-f option) <br>
- infers the N most-frequent haplotypes (option: -p N) for each interval <br>
  - supports homopolymer compression (option: -c), where reads are first compressed before the haplotypes are estimated <br>

<br>

### Quick start <br>
Here's an example:

python3 pycision.py  -q 4 -p 2 ./PrecisionID_mtDNA_WG_targets.bed  10_S17_L001.sorted.bam

This will use a minimum mapping quality of 4, and return the top 2 most frequent haplotypes for each amplicon in the above bam file. The bam file should be sorted and indexed prior to use.
<br>
A separate bam file will be made that has been trimmed to just the reads that span the entirety of the amplicons described in the bed file.
<br>


### Code requirements

Pycision is written for Python (version 3.*) and it relies on the pysam library. The pysam module is avaiable here <br>
https://pypi.python.org/pypi/pysam
<br>
I honestly have my doubts as to whether it can work on a windows system (perhaps though?), but on a typical unix system it should install relatively easily. <br>
Also, a bed file (PrecisionID_mtDNA_WG_targets.bed) is required. This file gives the locations of the regions targetted (sans the primers, which are treated as an unknown). As per the paper, I would modify locus mt_31 to end at position 3106 instead of 3107 to better handle the N in the reference sequence.
<br>
The file: PrecisionID_mtDNA_WG_targets.modified.bed has the modified coordinate for locus mt_31 (3106), and I would recommend using this bed in future projects. <br>
Run: <br>
python3 pycision.py <br>
To get a description of how to use this program...
<br><br>


