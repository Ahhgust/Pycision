#!/usr/bin/env python3

import os
import sys
import pysam
import csv
import re
import copy

# global. from the sam file specificatio. bam codes that consume the reference sequence
# taken from: https://samtools.github.io/hts-specs/SAMv1.pdf
# c[i] iff the bedop consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True]
# ditto for consuming the query
consumesQuery = [True, True, False, False, True, False, False, True]


def die(message=''):
    sys.stderr.write( message + os.linesep + "Oops! Correct usage:" + os.linesep + sys.argv[0] + " bedFile bamFile1 (...)" + os.linesep )
    sys.exit(1)


# returns the 1-based indexing start/stop posititions of the read in the reference sequence
def getGenomeStartStop(read):

# 1-based indexing
    currPos = genomeStart = read.pos + 1

    for (cigarType, cigarLength) in read.cigartuples:
        # M, D, or N (or = or X) (all consume reference bases)
        if (consumesReference[cigarType]):
            currPos += cigarLength

    return(genomeStart, currPos)

# performs soft clipping:
# changes the cigar string in READ to clip out regions outside of the bedRec (1-based indexed)
def makeSoftclippedCigar(read, genomePositions, bedRec):

    # number of bases in REFERENCE coordinates to the right to remove
    diffRight = genomePositions[1] - bedRec[2] 

# the old cigar string (converted from tuples to lists)
    ciggy = [list(elem) for elem in read.cigar ]

# should never happen...
    if diffRight <0:
        die("Unexpected trim: " + read + genomePositions + bedRec)

# trim to the right!    
    if diffRight > 0:
        numQConsumed = 0 # number of bases in the query (read) consumed

# pops off those items in the cigar that are (in there entirity) being converted into a soft-clop
        while ciggy:
            (t, l) = ciggy[-1] # type and length from the cigar
            
            if (consumesReference[t]):
                if l >= diffRight:
                    break
            
            ciggy.pop(-1)
            if consumesQuery[t]:
                numQConsumed += l

            if consumesReference[t]:
                diffRight -= l

            
        if (not ciggy):
            die("Unexpected cigar in read "  + read)


#only break out of the loop iff consumesReference[t] is true
        if consumesQuery[t]: # both are true
            numQConsumed += diffRight # so reduce this cigarop accordingly
        
        # shrink the (now) last op to the position of the reference sequence marked in the bed record
        ciggy[-1][1] -= diffRight

        # and if query bases were lost, soft clip them
        if (numQConsumed > 0):
            ciggy.append( list([4, numQConsumed]))
 


# and to the left!!
    diffLeft = bedRec[1] - genomePositions[0]
    
    if diffLeft > 0:
        numQConsumed = 0 # number of bases in the query (read) consumed

        while ciggy:
            (t, l) = ciggy[0] # type and length from the cigar
            
            if (consumesReference[t]):
                if l >= diffLeft:
                    break
            
            ciggy.pop(0)
            if consumesQuery[t]:
                numQConsumed += l

            if consumesReference[t]:
                diffLeft -= l

# no ciggy left. that's not good!
        if not ciggy:
            die("Unexpected cigar in read " + read)

# remember consumesRef is also true (by construction)
        if consumesQuery[t]:
            numQConsumed += diffLeft # residual differences added

        ciggy[0][1] -= diffLeft
        if numQConsumed > 0:
            ciggy.insert(0, list( [4, numQConsumed] ) )


#    print(ciggy, read)

# and convert things back to a tuple
    return( [tuple(elem) for elem in ciggy ] )




def softclipBam(bamFile, bedRecs):
    inBam = pysam.AlignmentFile(bamFile, "rb")
    outFile = re.sub('\.bam$', '.softClipped.bam', bamFile)
    outFileSorted = re.sub('\.bam$', '.softClipped.sorted.bam', bamFile)
    outBam = pysam.AlignmentFile(outFile, "wb", template=inBam)

    limits = [bedRecs[0][0], bedRecs[0][1], bedRecs[-1][2]]

    prevPos = -1
    currIndex = 0
    bedLen = len(bedRecs)
    # just look at reads in the mito. (the coordinates used/ mito padding used are immaterial with this approach)
    for read in inBam.fetch(limits[0], limits[1], limits[2]):
        if (read.is_unmapped):
            continue

        positions = getGenomeStartStop(read)

# TODO: adjust currIndex, and stop coordinates (assuming bam is sorted)        
        for i in range(currIndex, bedLen):
            rec = bedRecs[i]
            if (positions[0] <= rec[1] and positions[1] >= rec[2]):
                newcig = makeSoftclippedCigar(read, positions, rec)
                newread = copy.copy(read)
                newread.cigartuples = newcig
                newread.pos = rec[1] - 1 # soft clipping in front of read means we need to adjust the start coordinate. -1 b/c 0-based indexing in bam
                outBam.write(newread)
                break

    inBam.close()
    outBam.close()
# sort and index the output files.
    pysam.sort("-o", outFileSorted, outFile)
    pysam.index(outFileSorted, outFileSorted + ".bai")


if __name__ == "__main__":
    args = sys.argv
    if (len(args) < 3):
        die("I need at least one bed file and one bam file!")

# lazy person's bed file processing
    dat = open(args[1], 'r')
    bedRecs = []
    reader = csv.reader(dat, delimiter='\t')
    for row in reader:
        if (row[0][0] != '#'):
            row[1] = int(row[1]) + 1 # change from 0-based half open coordinates to 1-based
            row[2] = int(row[2])
            bedRecs.append(row)

# iterate over the bams, and soft clip them to just the regions
    for f in args[2:]:
        softclipBam(f, bedRecs)
        


