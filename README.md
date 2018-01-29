# pycision (V 0.001)

### Code requirements

Pycision is written for Python (version 3.*) and it relies on the pysam library. The pysam module is avaiable here <br>
https://pypi.python.org/pypi/pysam
<br>
I honestly have my doubts as to whether it can work on a windows system (perhaps though?), but on a typical unix system it should install relatively easily. <br>
Also, a bed file (PrecisionID_mtDNA_WG_targets.bed) is required. This file gives the locations of the regions targetted (sans the primers, which are treated as an unknown). As per the paper, I would modify locus mt_31 to end at position 3106 instead of 3107 to better handle the N in the reference sequence.
<br>
Run: <br>
python3 pycision.py <br>
To get a description of how to use this program...

