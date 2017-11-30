#!/usr/bin/env python3
# Written by August Woerner.
# This program takes in a bed file (arg1), and a bunch of sorted bams
# and it takes the intervals in the bed, and outputs a bunch of sorted, indexed bams
# that only contain reads that span at least one of the intervals in the bed.
# the reads are soft-clipped s.t. they correspond to *just* the positions specied.
# In practice, you give it a bed file of places you want to look at that have been PCR'd
# and an alignment that may include some/all of the primers
# and this reduces everything to just the regions sequenced by the primers w/o including the primer
# sequence itself.


import os
import sys
import pysam
import csv
import re
import copy
import argparse


# global. from the sam file specificatio. bam codes that consume the reference sequence
# taken from: https://samtools.github.io/hts-specs/SAMv1.pdf
# c[i] iff the bedop consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True]
# ditto for consuming the query
consumesQuery = [True, True, False, False, True, False, False, True]

VERSION = 0.0001

def die(message=''):
    sys.stderr.write( message + os.linesep + "Version number: " + str(VERSION) + os.linesep + "Oops! Correct usage:" + os.linesep + sys.argv[0] + " bedFile bamFile1 (...)" + os.linesep )
    sys.exit(1)


# returns the 1-based indexing start/stop posititions of the read in the reference sequence
def getGenomeStartStop(read):

# 1-based indexing
    currPos = genomeStart = read.pos + 1

    for (cigarType, cigarLength) in read.cigartuples:
        # M, D, or N (or = or X) (all consume reference bases)
        if (consumesReference[cigarType]):
            currPos += cigarLength

    return(genomeStart, currPos-1)

# performs soft clipping:
# changes the cigar string in READ to clip out regions outside of the bedRec (1-based indexed)
# of  trimLeft and trimRight, at least one must be true!
def makeSoftclippedCigar(read, genomePositions, bedRec, trimLeft, trimRight):

    # number of bases in REFERENCE coordinates to the right to remove
    diffRight = genomePositions[1] - bedRec[2] 

# the old cigar string (converted from tuples to lists)
    ciggy = [list(elem) for elem in read.cigar ]

# should never happen...
    if diffRight <0 and trimRight:
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
        

        if diffRight < ciggy[-1][1]:
        # shrink the (now) last op to the position of the reference sequence marked in the bed record
            ciggy[-1][1] -= diffRight
        else: # avoid creating 0-length cigar ops
            ciggy.pop(-1)

        # and if query bases were lost, soft clip them
        if (numQConsumed > 0):
            ciggy.append( list([4, numQConsumed]))
 


# and to the left!!
    diffLeft = bedRec[1] - genomePositions[0]
    if diffLeft <0 and trimLeft:
        die("Unexpected trim: " + read + genomePositions + bedRec)

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

        if diffLeft < ciggy[0][1]:
            ciggy[0][1] -= diffLeft
        else:
            ciggy.pop(0)

        if numQConsumed > 0:
            ciggy.insert(0, list( [4, numQConsumed] ) )


#    print(ciggy, read)

# and convert things back to a tuple
    return( [tuple(elem) for elem in ciggy ] )




def softclipBam(bamFile, bedRecs, halfway):
    inBam = pysam.AlignmentFile(bamFile, "rb")
    outFile = re.sub('\.bam$', '.softClipped.bam', bamFile)
    outFileSorted = re.sub('\.bam$', '.softClipped.sorted.bam', bamFile)
    outBam = pysam.AlignmentFile(outFile, "wb", template=inBam)

    limits = [bedRecs[0][0], bedRecs[0][1], bedRecs[-1][2]]

    prevPos = -1
    currIndex = 0
    bedLen = len(bedRecs)
    previousStart = -1
    # just look at reads in the mito. (the coordinates used/ mito padding used are immaterial with this approach)
    for read in inBam.fetch(limits[0], limits[1], limits[2]):
        if (read.is_unmapped):
            continue

        if (read.pos < previousStart):
            die("File: " + bamFile + " does not appear to be sorted!" + read.pos + " vs " + previousStart)

        previousStart = read.pos
        positions = getGenomeStartStop(read)

        while (currIndex < bedLen and read.pos > bedRecs[currIndex][2]):
            currIndex += 1

        for i in range(currIndex, bedLen):
            rec = bedRecs[i]
            if (positions[0] <= rec[1] and positions[1] >= rec[2]):
                newcig = makeSoftclippedCigar(read, positions, rec, True, True)
                newread = copy.copy(read)
                newread.cigartuples = newcig
                newread.pos = rec[1] - 1 # soft clipping in front of read means we need to adjust the start coordinate. -1 b/c 0-based indexing in bam
                outBam.write(newread)
#                print("Both")
                break
            elif (halfway and positions[0] <= rec[1] and positions[1] >= rec[1]/2.0 + rec[2]/2.0): # to the left, to the left...
                newcig = makeSoftclippedCigar(read, positions, rec, True, False)
                newread = copy.copy(read)
                newread.cigartuples = newcig
                newread.pos = rec[1] - 1 # soft clipping in front of read means we need to adjust the start coordinate. -1 b/c 0-based indexing in bam
                outBam.write(newread)
 #               print("Left")
                break
            elif (halfway and positions[1] >= rec[2] and positions[0] <= rec[1]/2.0 + rec[2]/2.0): # to the right, to the right... everybody sing
                newcig = makeSoftclippedCigar(read, positions, rec, False, True)
                newread = copy.copy(read)
                newread.cigartuples = newcig
                #newread.pos = rec[1] - 1 # the right-hand position wasn't trimmed, so no need to adjust.
                outBam.write(newread)
  #              print("Right")
                break
            elif positions[1] < rec[2]:
                break

   
    inBam.close()
    outBam.close()
    # sort and index the output files.
    pysam.sort("-o", outFileSorted, outFile)
    pysam.index(outFileSorted, outFileSorted + ".bai")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Let's trim some bams!")

    halfway=False
    parser.add_argument('-f', '--fifty_percent', dest='halfway', help="Keeps reads that span to halfway through a region", action='store_true')

    results = parser.parse_known_args(sys.argv[1:])[0]
    halfway = results.halfway

    if halfway:
        args = sys.argv[2:]
    else:
        args = sys.argv[1:]


    if (len(args) < 2):
        parser.print_help()
        die("I need at least one bed file and one bam file!")

# lazy person's bed file processing
    dat = open(args[0], 'r')
    bedRecs = []
    reader = csv.reader(dat, delimiter='\t')
    chrom =''

    for row in reader:
        if (row[0][0] != '#'):
            row[1] = int(row[1]) + 1 # change from 0-based half open coordinates to 1-based
            row[2] = int(row[2])
            if chrom == '':
                chrom = row[0]
            elif chrom != row[0]:
                parser.print_help()
                die("At least two separate chromosomes detected\n" + chrom + "\nAnd\n" + row[0] + "\nThat's not good!!\n")

            bedRecs.append(row)

# iterate over the bams, and soft clip them to just the regions
    for f in args[1:]:
        print(f)
        softclipBam(f, bedRecs, halfway)
        


