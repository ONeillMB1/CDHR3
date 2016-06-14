#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt

#######################################################################
##This script reads in a ped file and divides it into multiple ped  
##files divided by population.
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
    parser.add_argument('-s', '--sample',
        help = 'file containing paths to bamqc results that are to be collated', 
        #type = is_dir,
        required = True)
    return parser.parse_args()

def pop_dict():
    d = {}
    p = {}
    with open(args.sample, 'r') as pops:
        for line in pops:
            if line[0] != "#":
                line = line.strip()
                samp,pop,superpop,gender = line.split('\t')
                d[samp] = [pop,superpop]
                if pop in p:
                    continue
                else:
                    p[pop] = []
                if superpop in p:
                    continue
                else:
                    p[superpop] = []
    return d,p

def split_ped(d,p):
    with open(args.ped, 'r') as ped:
        for line in ped:
            line = line.strip()
            samp = line.split()[0]
            p[d[samp][0]].append(line)
            p[d[samp][1]].append(line)
    prefix = args.ped.split(".")[0]
    for pop in p:
        print("There are {0} subjects in population {1}".format(len(p[pop]), pop))
        if len(p[pop]) > 1:
            outfile = pop + '_' + prefix + '.ped'
            print("Writing info to {0}".format(outfile))
            with open(outfile, 'w') as out:
                for i in p[pop]:
                    out.write(i + '\n')
            

    

args = get_arguments()
d,p = pop_dict()
split_ped(d,p)
