
import re;
import os;
import sys;
import string;
import math;

import numpy;
import peakutils;
import argparse;
from argparse import RawTextHelpFormatter

import logging

import findTrinucleotideRepeats

from myheader import *

def non_negative(i, mstr):
	if i<0: 
		print (("Error %d could not be negative(%d)" % (mstr, i)))
		sys.exit(1);



parser = argparse.ArgumentParser(description='''Find trinucleotide or hexanucleotide 7+ repeats for (a) gene(s) of interests or for all genes in the test genome of the BAM file (for hg38!!!!!!!!).''', epilog="For example, \n \
\tpython %(prog)s freeze4-all-merge.sort.bam; \n \
\tpython %(prog)s freeze4-all-merge.sort.bam -repeatgene HTT; \n \
\tpython %(prog)s freeze4-all-merge.sort.bam -repeatgene ATXN3 --align 1 --updown 18 --extend 0; \n \
OR \tpython %(prog)s freeze4-all-merge.sort.bam -repeatgene all; \n \
Final results was stored in TrinRepDis.log", formatter_class=RawTextHelpFormatter);
parser.add_argument("bamfile", help="A BAM file storing all alignments");
parser.add_argument("-hg", default='hg38', help="The reference genome is used. Currently, only hg38 is supported");
parser.add_argument("--align", type=int, default=1, help="Is unsymmetrical alignment used for error correction. 1: yes, 0: no");
parser.add_argument("--updown", type=int, default=18, help="Is upstream/downstream used for repeat length inference. non-0: yes, 0: no");
parser.add_argument("--extend", type=int, default=0, help="Is upstream/downstream extended as repeat region. non-0: yes, 0: no");
parser.add_argument("-repeatgene", help="A gene name which you want to analyze, such as HTT for huntington's disease. \
                                        'all': all known genes associated with trinucleotide repeat disorders will be analyzed;");

#LOG_FILENAME = 'TrinRepDis.log'

args = parser.parse_args();

logfolder = 'log/'

if not os.path.isdir(logfolder):
	os.system('mkdir '+logfolder)

isAlign, isupdown, isExtend = args.align, args.updown, args.extend
non_negative(isAlign, 'isAlign');
non_negative(isupdown, 'isupdown');
non_negative(isExtend, 'isExtend');

simulation_file_id = ""
analysis_file_id = ""

if isAlign>0: analysis_file_id += '_align'
else: analysis_file_id += '_non_align'

if isupdown>0: analysis_file_id += '_updn'+str(isupdown)
else: analysis_file_id += '_non_updn'

if isExtend>0: analysis_file_id += '_ext'+str(isExtend)
else: analysis_file_id += '_non_ext'

analysis_file_id += ("for_%s" % (sys.argv[1]).replace('/', '_'))

simulation_file_id = ""
unique_file_id = simulation_file_id+analysis_file_id

if args.repeatgene==None or args.repeatgene=='all':
        repeatgene = 'all'
else: repeatgene = args.repeatgene

filename = logfolder + 'TrinRepDis_' + repeatgene + unique_file_id + '.log'

LOG_FILENAME = filename
logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")


logging.info('Input BAM file is '+args.bamfile);
gLoc = findTrinucleotideRepeats.getDiseseGeneInRefGenomeLocation(args.hg)
summary = []
#print args.repeatgene
if args.repeatgene==None or args.repeatgene=='all':
	summary = findTrinucleotideRepeats.getRepeat(gLoc, args.bamfile, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id)
else:
	summary.append(findTrinucleotideRepeats.getRepeatForKnownGene(gLoc, args.repeatgene, args.bamfile, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id))

for curgrep in summary:
	pstr = 'NONE'
	if len(curgrep[2])>0: pstr = str(curgrep[2][0])
	if len(curgrep[2])>1: pstr += ','+ str(curgrep[2][1])
	logging.critical('\t Gene name='+curgrep[0]+'; ref_repeat='+('%.0f' % (curgrep[1]))+('(allreads=%d); ' % (curgrep[4]))+'repeat in test genome='+ pstr + ' '+str(curgrep[5]))
	logging.critical('\t#Format: Repeat_time=numberOfAlignedRead:'+curgrep[3])
	

