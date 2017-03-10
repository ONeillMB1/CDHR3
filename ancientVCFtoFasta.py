#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt

#######################################################################
##This script reads in a vcf file and converts it to a fasta accounting  
##for sequence gaps with N.
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
    parser.add_argument('-v', '--vcf',
        help = 'vcf file', 
        type = is_file,
        required = True)
    return parser.parse_args()

def split_ind():
    d = {}
    with open(args.vcf, 'r') as vcf:
        for line in vcf:
            if line[0] != "#":
                line = line.strip().split('\t')
                pos = int(line[1])
                ref = line[3]
                alt = line[4]
                gt = line[9].split(":")[0]
                GT = gt.replace("0", ref)
                GT = GT.replace("1", alt)
                GT = GT.split("/")
                d[pos] = GT
    dd = {}
    for key in d:
        if len(set(d[key])) > 1:
            print(str(key) + '\t' + ",".join(d[key]))
            dd[key] = {'N'}
        else:
            dd[key] = set(d[key])
    return(dd)

def write_haplotype(dd):
    sequence = []
    for i in range(105657078, 105680023):
        if i in dd:
            allele = next(iter(dd[i]))
            sequence.append(allele)
        else:
            sequence.append("-")
    outfilename = args.vcf.split(".")[0] + ".fasta"
    with open(outfilename, 'w') as outfile:
        outfile.write(">" + args.vcf.split("_")[0] + '\n')
        outfile.write("".join(sequence) + '\n')
    return(sequence)


args = get_arguments()
dd = split_ind()
sequence = write_haplotype(dd)
