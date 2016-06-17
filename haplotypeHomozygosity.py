#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt

#######################################################################
##This script reads in a ped file and converts it to a snp alignment  
##(fasta)
#######################################################################


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_arguments(): 
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Collate summary info from multiple bamqc runs")
    parser.add_argument('-p', '--ped',
        help = 'plink ped file', 
        type = is_file,
        required = True)
    parser.add_argument('-m', '--map',
        help = 'plink map file', 
        type = is_file,
        required = True)
    return parser.parse_args()

def find_index(lst):
    return [i for i, x in enumerate(lst) if x == 105658451]

def get_pos():
    pos = []
    with open(args.map, 'r') as mapFile:
        for line in mapFile:
            if line[0] != "#":
                line = line.strip().split()
                position = int(line[3])
                pos.append(position)
    rev_pos = pos[::-1]
    return(pos, rev_pos)

def hap_hom(pos, rev_pos):
    outfileName = args.ped.split(".")[0] + "_haplotypeHomozygoxity.txt"
    with open(args.ped, 'r') as pedFile, open(outfileName, 'w') as outFile:
        for line in pedFile:
            if line[0] != "#":
                line = line.strip().split()
                sample = line[0]
                gt = zip(line[6::2], line[7::2])
                rev_gt = gt[::-1]
                idx = find_index(pos)[0]
                rev_idx = find_index(rev_pos)[0]
                allele = "".join(gt[idx])
                if gt[idx][0] == gt[idx][1]:
                    start = "-Inf"
                    stop = "Inf"
                    for i, pair in enumerate(gt[(idx+1):]):
                        if pair[0] != pair[1]:
                            stopIdx = i + idx + 1
                            stop = pos[stopIdx]
                            break
                        else:
                            continue
                    for i, pair in enumerate(rev_gt[(rev_idx+1):]):
                        if pair[0] != pair[1]:
                            startIdx = i + rev_idx + 1
                            start = rev_pos[startIdx]
                            break
                        else:
                            continue
                    outFile.write(sample + '\t' + allele + '\t' + str(start) + '\t' + str(stop) + '\n')
                else:
                    print("Individual {0} is heterozygous at rs6967730".format(sample))
                    outFile.write(sample + '\t' + allele + '\t' + "NA\tNA\n")
                    


args = get_arguments()
pos, rev_pos = get_pos()
hap_hom(pos, rev_pos)
