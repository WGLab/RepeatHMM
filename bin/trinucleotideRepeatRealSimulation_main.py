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

from scripts import trinucleotideRepeatRealSimulation
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
\tpython %(prog)s -repeatgene ATXN3 --align 1 --updown 18 --extend 0 --coverage 100 --randTimes 100 -UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1; \n \
OR \tpython %(prog)s -repeatgene all;\n \
Final results was stored in logsim/TrinRepSim*.log", formatter_class=RawTextHelpFormatter);
parser.add_argument("-hg", default='hg38', help="The reference genome is used. Currently, only hg38 is supported");
parser.add_argument("--align", type=int, default=1, help="Is unsymmetrical alignment used for error correction. 1: yes(Default), 0: no");
parser.add_argument("--updown", type=int, default=18, help="Is upstream/downstream used for repeat length inference. non-0: yes(Default: 18), 0: no");
parser.add_argument("--extend", type=int, default=0, help="Is upstream/downstream extended as repeat region. non-0: yes, 0: no(Default)");

parser.add_argument("-repeatgene", default=None, help="A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. \
                                        'all': all known genes associated with trinucleotide repeat disorders will be analyzed;");

parser.add_argument("--insert_rate", type=float, default=0.12, help="Insert error rate. Default: 0.12");
parser.add_argument("--del_rate", type=float, default=0.02, help="Deletion error rate. Default: 0.02");
parser.add_argument("--sub_rate", type=float, default=0.02, help="Substitution error rate. Default: 0.02");
parser.add_argument("--coverage", type=int, default=300, help="The number of reads produced after mutations. Default: 300");

parser.add_argument("--randTimes", type=int, default=100, help="The number of simulation times. Default: 100");
parser.add_argument("--UserDefinedGene", default=UserDefinedGenedefault, help="The gene information defined by users. If this option is given, the default gene information will be revised. Default: ///////");
parser.add_argument("--UserDefinedGeneName", default="Pcr1", help="The name for storing results. Default: Pcr1");
parser.add_argument("--UnsymAlign", type=int, default=1, help="Whether unsymmetrical alignment (Default: 1) rather than bwa mem (0) would be used");


args = parser.parse_args();


if args.repeatgene==None:
        parser.print_help()
        sys.exit(140)


logfolder = 'logsim/'

if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)

isAlign, isupdown, isExtend = args.align, args.updown, args.extend
non_negative(isAlign, 'isAlign');
non_negative(isupdown, 'isupdown');
non_negative(isExtend, 'isExtend');

isUnsymAlign = args.UnsymAlign

simulation_file_id = ""
analysis_file_id = ""

if isUnsymAlign: analysis_file_id += "_unsymalign"

if isAlign>0: analysis_file_id += '_align'
else: analysis_file_id += '_non_align'

if isupdown>0: analysis_file_id += '_updn'+str(isupdown)
else: analysis_file_id += '_non_updn'

if isExtend>0: analysis_file_id += '_ext'+str(isExtend)
else: analysis_file_id += '_non_ext'

insert_rate=args.insert_rate;
del_rate=args.del_rate
sub_rate=args.sub_rate
coverage = args.coverage
randTimes = args.randTimes

if args.repeatgene==None:#  or args.repeatgene=='all':
	#repeatgene = 'all'
	repeatgene = args.UserDefinedGeneName
else: repeatgene = args.repeatgene

#simulation_file_id += ('_ins%.2f_del%.2f_sub%.2f_cov%d_times%d' % (insert_rate, del_rate, sub_rate, coverage, randTimes))

simulation_file_id += ('_ins%.2f_del%.2f_sub%.2f_cov%d' % (insert_rate, del_rate, sub_rate, coverage))
analysis_file_id += ('_times%d' % randTimes)

if not args.UserDefinedGene == UserDefinedGenedefault:
	simulation_file_id = ("_def%s" % args.UserDefinedGeneName) + simulation_file_id

unique_file_id =  simulation_file_id+analysis_file_id

filename = logfolder + 'simTrinRepDis_' + repeatgene + unique_file_id + '.log'

#print filename, '\n', repeatgene,'\n',  unique_file_id, '\n', simulation_file_id, '\n', analysis_file_id, '\n', args.UserDefinedGene, '\n', args.UserDefinedGeneName

LOG_FILENAME = filename
logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

gLoc = findTrinucleotideRepeats.getDiseseGeneInRefGenomeLocation(args.hg)

random.seed(7);

if not args.UserDefinedGene == UserDefinedGenedefault:
	res = trinucleotideRepeatRealSimulation.getSimforKnownGeneWithPartialRev(gLoc, isUnsymAlign, unique_file_id, simulation_file_id, analysis_file_id, repeatgene, args.UserDefinedGene, insert_rate, del_rate, sub_rate, coverage, isAlign, isupdown, isExtend, randTimes)
else:
	if args.repeatgene=='all':
		res = trinucleotideRepeatRealSimulation.getSim(gLoc, isUnsymAlign, unique_file_id, simulation_file_id, analysis_file_id, insert_rate, del_rate, sub_rate, coverage, isAlign, isupdown, isExtend, randTimes);
	else:
		res = trinucleotideRepeatRealSimulation.getSimforKnownGene(gLoc, isUnsymAlign, unique_file_id, simulation_file_id, analysis_file_id, args.repeatgene, insert_rate, del_rate, sub_rate, coverage, isAlign, isupdown, isExtend, randTimes)

for r1 in res:
	if not isUnsymAlign:
		logging.info('%s, %s<<<True<<<>>>Pred>>>%s, %s \t\t#wrongAlign %s' % (str(r1[0][0]), str(r1[0][1]), str(r1[1][0]), str(r1[1][1]), str(r1[2])));
	else: logging.info('%s, %s<<<True<<<>>>Pred>>>%s, %s ' % (str(r1[0][0]), str(r1[0][1]), str(r1[1][0]), str(r1[1][1])));
	#print r1[0][0], r1[0][1], '<<<True<<<>>>Pred>>>', r1[1][0], r1[1][1], '\t\t#wrongAlign', r1[2];

simresfolder = 'sim_res/'
if not os.path.isdir(simresfolder):
	os.system('mkdir '+simresfolder)

curfilename = simresfolder + repeatgene + unique_file_id + '.txt'
curres = []
for r1 in res:
	curstr = ('%s %s %s' % (str(r1[0][0]), str(r1[0][1]), str(r1[1][0])))
	if len(r1[1])>1 and (not string.strip(str(r1[1][1]))==''): curstr += (' %s' % str(r1[1][1]))
	else: curstr += (' %s' % str(r1[1][0]))
	curres.append(curstr)
findTrinucleotideRepeats.myWriteTxtFile(curres, curfilename);

