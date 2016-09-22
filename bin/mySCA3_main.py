import re;
import os;
import sys;
import string;
import math;
import random;

import numpy;
import peakutils;
import argparse;

import logging

from argparse import RawTextHelpFormatter

from scripts import mySCA3
from scripts import findTrinucleotideRepeats

from scripts.myheader import *


def non_negative(i, mstr):
        if i<0:
                print (("Error %d could not be negative(%d)" % (mstr, i)))
                sys.exit(1);



#LOG_FILENAME = 'TrinRepSim.log'
#logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

UserDefinedGenedefault = "///////"

parser = argparse.ArgumentParser(description='''Simulate Trinucleotide repeat for (a) known gene(s) of interests or for all genes (for hg38!!!!!!!!).''', epilog="For example, \n \
\tpython %(prog)s ; \n \
\tpython %(prog)s -repeatgene HTT; \n \
\tpython %(prog)s -repeatgene ATXN3 --align 1 --updown 18 --extend 0; \n \
\tpython %(prog)s -repeatgene ATXN3 --align 1 --updown 18 --extend 0 -UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1; \n \
OR \tpython %(prog)s -repeatgene all;\n \
Final results was stored in logsca3/TrinRepSca3*.log", formatter_class=RawTextHelpFormatter);
parser.add_argument("-hg", default='hg38', help="The reference genome is used. Currently, only hg38 is supported");
parser.add_argument("--align", type=int, default=1, help="Is unsymmetrical alignment used for error correction. 1: yes(default), 0: no");
parser.add_argument("--updown", type=int, default=18, help="Is upstream/downstream used for repeat length inference. non-0: yes(default: 18), 0: no");
parser.add_argument("--extend", type=int, default=0, help="Is upstream/downstream extended as repeat region. non-0: yes, 0: no(default)");

parser.add_argument("-repeatgene", default=None, help="A gene name which you want to analyze, such as HTT for huntington's disease. Default: None");

parser.add_argument("--UserDefinedGene", default=UserDefinedGenedefault, help="The gene information defined by users. If this option is given, the default gene information will be revised. Default: ///////");
parser.add_argument("--UserDefinedGeneName", default="sca3_Pcr1", help="The name for storing results. Default: sca3_Pcr1");
parser.add_argument("--UnsymAlign", type=int, default=1, help="Whether unsymmetrical alignment (1) rather than bwa mem (0) would be used. Default: 1");
parser.add_argument("--fastq", help="The file name for fasta sequences");

args = parser.parse_args();

if args.repeatgene==None:
	parser.print_help()
	sys.exit(140)

logfolder = 'logsca3/'

if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)

isAlign, isupdown, isExtend = args.align, args.updown, args.extend
non_negative(isAlign, 'isAlign');
non_negative(isupdown, 'isupdown');
non_negative(isExtend, 'isExtend');

isUnsymAlign = args.UnsymAlign

analysis_file_id = ""

if isUnsymAlign: analysis_file_id += "_unsymalign"

if isAlign>0: analysis_file_id += '_align'
else: analysis_file_id += '_non_align'

if isupdown>0: analysis_file_id += '_updn'+str(isupdown)
else: analysis_file_id += '_non_updn'

if isExtend>0: analysis_file_id += '_ext'+str(isExtend)
else: analysis_file_id += '_non_ext'

if args.repeatgene==None:#
	repeatgene = args.UserDefinedGeneName
else: repeatgene = args.repeatgene

unique_file_id =  args.UserDefinedGeneName+analysis_file_id

fastafile=args.fastq;

filename = logfolder + 'sca3TrinRepDis_' + repeatgene + unique_file_id + '.log'

LOG_FILENAME = filename
logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

gLoc = findTrinucleotideRepeats.getDiseseGeneInRefGenomeLocation(args.hg)

res = mySCA3.getSCA3forKnownGeneWithPartialRev(gLoc, fastafile, isUnsymAlign, unique_file_id, analysis_file_id, repeatgene, args.UserDefinedGene, isAlign, isupdown, isExtend)


for r1 in res:
	print r1
	#print r1[:2], r1[2];


