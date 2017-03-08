#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt
import glob
import subprocess
import shlex

####################################################################
##This script reads in the header of a vcf file to extract individual
##identifiers and then looks at the SNP and returns a list of those
##individuals that are heterozygous.
####################################################################


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname
        
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
        description="Merge individual sample VCFs into a project VCF")
    parser.add_argument('-v', '--vcf', required=True,
        help = 'vcf file',
        type = is_file)
    parser.add_argument('-o', '--out', required=True,
        help = 'out file')
    return parser.parse_args()

args = get_arguments()

def find_hets():
    with open(args.vcf, 'r') as vcf, open(args.out, 'w') as outfile:
        for line in vcf:
            line = line.strip()
            if line[0] == "#":
                if line[0:6] == "#CHROM":
                    line = line.split('\t')
                    ids = line[9:]
                else:
                    continue
            else:
                hets = []
                line = line.split('\t')
                snpID = "_".join(line[0:3])
                print snpID
                genotypes = line[9:]
                for i, gen in enumerate(genotypes):
                    if gen == "0|1" or gen == "1|0":
                        hets.append(i)
                outfile.write(snpID + '\n')
                for x in hets:
                    outfile.write(ids[x] + '\n')
    return ids, genotypes

ids, genotypes = find_hets()
