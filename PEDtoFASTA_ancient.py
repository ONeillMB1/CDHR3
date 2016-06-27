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
        help = 'file containing paths to bamqc results that are to be collated', 
        type = is_file,
        required = True)
    return parser.parse_args()

def split_ind():
    outfile = args.ped.split(".")[0] + ".fasta"
    with open(args.ped, 'r') as ped, open(outfile, 'w') as outfile:
        for line in ped:
            if line[0] != "#":
                line = line.strip().split()
                sample = line[0]
                #prehap1 = line[6::2]
                hap1 = [x.replace('0', '?') for x in line[6::2]]
                #prehap2 = line[7::2]
                hap2 = [x.replace('0', '?') for x in line[7::2]]
                if hap1 == hap2:
                    missing = float(hap1.count("?"))/len(hap1)
                    if missing < 0.25:
                        outfile.write(">" + sample + "-" + hap1[0] + "\n")
                        outfile.write("".join(hap1) + '\n')
                else:
                    outfile.write(">" + sample + "_h1-" + hap1[0] + "\n")
                    outfile.write("".join(hap1) + '\n')
                    outfile.write(">" + sample + "_h2-" + hap2[0] + "\n")
                    outfile.write("".join(hap2) + '\n')
args = get_arguments()
split_ind()
