#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC,Gapped
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
        description="Fill in unknown sites for ancient sample")
    parser.add_argument('-r', '--refMap',
        help = 'map file of dataset to be combined with', 
        type = is_file,
        required = True)
    parser.add_argument('-a', '--ancMap',
        help = 'map file of ancient data', 
        type = is_file,
        required = True)
    parser.add_argument('-f', '--ancFasta',
        help = 'fasta corresponding to ancient map file', 
        type = is_file,
        required = True)
    return parser.parse_args()

def find_idx():
    refSNPs = []
    with open(args.refMap, 'r') as ref:
        for line in ref:
            if line[0] != "#":
                line = line.strip().split()
                chrom = line[0]
                pos = line[3]
                snpID = chrom + '_' + pos
                refSNPs.append(snpID)
    print len(refSNPs)
    ancSNPs = []
    with open(args.ancMap, 'r') as anc:
        for line in anc:
            if line[0] != "#":
                line = line.strip().split()
                chrom = line[0]
                pos = line[3]
                snpID = chrom + '_' + pos
                ancSNPs.append(snpID)
    print len(ancSNPs)
    indices = []
    for i in ancSNPs:
        try:
            idx = refSNPs.index(i)
        except ValueError:
            idx = "NA"
        indices.append(idx)
    return refSNPs, indices

def write_fasta(refSNPs, indices):
    outfileName = args.ancFasta.split(".")[0] + "_toBeCombined.fasta"
    records = []
    with open(outfileName, 'w') as outfile: 
        for seq_record in SeqIO.parse(args.ancFasta, "fasta"):
            if len(seq_record.seq) == len(indices):
                newRec = ['?']*len(refSNPs)
                for i, idx in enumerate(indices):
                    if idx != "NA":
                        newRec[idx] = seq_record.seq[i]
                    else:
                        continue
                outfile.write('>' + seq_record.id + '\n')
                outfile.write("".join(newRec) + '\n')

args = get_arguments()
refSNPs, indices = find_idx()
write_fasta(refSNPs, indices)
