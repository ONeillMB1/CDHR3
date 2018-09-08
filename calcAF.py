#!/usr/bin/env python

import argparse
import sys
from Bio import SeqIO
import os
import getopt

#######################################################################
##This script reads in a vcf file and a file listing the individuals 
##and there respective populations and calculates the allele frequency
##at each site in the vcf for each population.
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
        description="calculate allele frequencies")
    parser.add_argument('-v', '--vcf',
        help = 'vcf file', 
        type = is_file,
        required = True)
    parser.add_argument('-i', '--ind',
        help = 'ind file', 
        type = is_file,
        required = True)
    parser.add_argument('-o', '--out',
        help = 'out file', 
        required = True)
    return parser.parse_args()

def get_pops():
    d = {}
    pops = set()
    with open(args.ind, 'r') as indFile:
        for line in indFile:
            if line[0] != "#":
                line = line.strip().split('\t')
                if line[7] != 'IGNORE':
                    d[line[0]] = line[7]
                    pops.add(line[7])
                else:
                    continue
    return d, pops

def calc_AF(d, pops):
    popOrd = sorted(pops)
    with open(args.vcf, 'r') as vcfFile, open(args.out, 'w') as outfile:
        for line in vcfFile:
            line = line.strip()
            if line[0] == "#":
                if 'CHROM' in line:
                    header = line.split('\t')
                    samples = header[9:]
                    outfile.write("\t".join(header[0:5]) + '\t' + "\t".join(popOrd) + '\t' + "total" + '\n')
                    #print(",".join(samples))
                else:
                    continue
            else:
                alleles = []
                AD = {}
                for pop in pops:
                    AD[pop] = []
                line = line.split('\t')
                genotypes = line[9:]
                for i, gt in enumerate(genotypes):
                    if samples[i] in d:
                        pop = (d[samples[i]])
                        gts = gt.split("/")
                        if gts[0] == gts[1]:
                            if str(gts[0]) == '0':
                                allele = line[3]
                                alleles.append(allele)
                                AD[pop].append(allele)
                            elif str(gts[0]) == '1':
                                allele = line[4]
                                alleles.append(allele)
                                AD[pop].append(allele)
                            else:
                                continue
                #print(" ".join(alleles))
                try:
                    #refFreq = alleles.count(line[3])/float(len(alleles))
                    altFreq = alleles.count(line[4])/float(len(alleles))
                except ZeroDivisionError:
                    #refFreq = "NA"
                    altFreq = "NA"
                #print("Reference allele ({0}) has frequency {1}".format(line[3],refFreq))
                #print("Alternative allele ({0}) has frequency {1}".format(line[4],altFreq))
                outfile.write("\t".join(line[0:5]) + '\t')
                for x in popOrd:
                    try:
                        #afR = AD[x].count(line[3])/float(len(AD[x]))
                        afA = AD[x].count(line[4])/float(len(AD[x]))
                    except ZeroDivisionError:
                        afA = "NA"
                    outfile.write(str(afA) + '\t')
                    #print(pop + '\t' + line[3] + '\t' + str(afR))
                    #print(pop + '\t' + line[4] + '\t' + str(afA))
                outfile.write(str(altFreq) + '\n')

    return AD

args = get_arguments()
d, pops = get_pops()
AD = calc_AF(d, pops)
