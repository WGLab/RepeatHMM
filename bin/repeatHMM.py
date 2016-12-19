
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

#from scripts import trinucleotideRepeatRealSimulation
#from scripts import findTrinucleotideRepeats
#from scripts.myheader import *

from scripts import trinucleotideRepeatRealSimulation
from scripts import findTrinucleotideRepeats
from scripts import mySCA3

from scripts.myheader import *


parser = argparse.ArgumentParser(description="Determine Trinucleotide repeat for (a) gene(s) of interests or for all genes.", epilog="For example, \n \
\tpython %(prog)s BAMinput: with a BAM file as input\n \
\tpython %(prog)s FASTQinput: with a FASTQ file as input \n \
\tpython %(prog)s Simulate: for simulation", formatter_class=RawTextHelpFormatter);


originalError = '!!!Error: !!!!!! \n'


def non_negative(i, mstr):
   if i<0: return (("\tError %d could not be negative(%d)\n" % (mstr, i)))
   else: return ''

def getCommonOptions(margs):
   errorStr = ""
  
   isRemInDel, isupdown, isExtend = margs.RemInDel, margs.updown, margs.extend
   errorStr += non_negative(isRemInDel, 'isRemInDel');
   errorStr += non_negative(isupdown, 'isupdown');
   errorStr += non_negative(isExtend, 'isExtend');

   SeqDepth = margs.SeqDepth;
   errorStr += non_negative(SeqDepth, 'SeqDepth');

   analysis_file_id_common = ''
   if isRemInDel>0: analysis_file_id_common += '_RemInDel'
   else: analysis_file_id_common += '_non_RemInDel'

   if isupdown>0: analysis_file_id_common += '_updn'+str(isupdown)
   else: analysis_file_id_common += '_non_updn'

   if isExtend>0: analysis_file_id_common += '_ext'+str(isExtend)
   else: analysis_file_id_common += '_non_ext'

   if not os.path.isfile(hg_reference_and_index+'/'+margs.hg+'.fa'):
      errorStr += '\tNo reference genome under the folder '+hg_reference_and_index+' with the filename '+margs.hg+'.fa\n'
   gLoc = findTrinucleotideRepeats.getDiseseGeneInRefGenomeLocation(margs.hg)
   
   specifiedGeneName = margs.UserDefinedGeneName
   knownGeneName = margs.repeatgene
   specifiedGeneInfo = margs.UserDefinedGene

   if UserDefinedGenedefault==specifiedGeneInfo and knownGeneName==None:
      errorStr += '\tNone gene information is given. \n'
   if (not specifiedGeneInfo==UserDefinedGenedefault) and knownGeneName=='all':
      errorStr += '\t"repeatgene" cannot be "all" when specifying "UserDefinedGene"\n'
 
   commonOptions = {}
   commonOptions['SeqDepth'] = SeqDepth
   commonOptions['isRemInDel'] = isRemInDel
   commonOptions['isupdown'] = isupdown
   commonOptions['isExtend'] = isExtend
   commonOptions['specifiedGeneName'] = specifiedGeneName
   commonOptions['knownGeneName'] = knownGeneName
   commonOptions['specifiedGeneInfo'] = specifiedGeneInfo
   if knownGeneName == None:
      if (not specifiedGeneName==None): 
         commonOptions['repeatgene'] = specifiedGeneName
      else: commonOptions['repeatgene'] = None;
   else: commonOptions['repeatgene'] = knownGeneName
   commonOptions['hg'] = margs.hg
   if margs.hgfile==None or margs.hgfile=='':
      commonOptions['hgfile'] = commonOptions['hg']+'.fa'
   else:
      commonOptions['hgfile'] = margs.hgfile
   commonOptions['gLoc'] = gLoc

   analysis_file_id_common += '_'+commonOptions['hg']

   return commonOptions, errorStr, analysis_file_id_common

def simulation(margs):
   specifiedOptions = {}
   simulation_file_id = ""
   analysis_file_id = ""

   errorStr = originalError
   commonOptions, errorStr_com, analysis_file_id_com = getCommonOptions(margs)
   errorStr += errorStr_com
  
   isUnsymAlign = margs.UnsymAlign; specifiedOptions['isUnsymAlign'] = isUnsymAlign
   if isUnsymAlign: analysis_file_id += "_unsymalign"
   analysis_file_id += analysis_file_id_com

   insert_rate=margs.insert_rate; errorStr += non_negative(insert_rate, 'insert_rate');
   del_rate=margs.del_rate;  errorStr += non_negative(del_rate, 'del_rate');
   sub_rate=margs.sub_rate;  errorStr += non_negative(sub_rate, 'sub_rate');
   coverage = margs.coverage; errorStr += non_negative(coverage, 'coverage');
   randTimes = margs.randTimes; errorStr += non_negative(randTimes, 'randTimes');

   specifiedOptions['insert_rate'] = insert_rate
   specifiedOptions['del_rate'] = del_rate
   specifiedOptions['sub_rate'] = sub_rate
   specifiedOptions['coverage'] = coverage
   specifiedOptions['randTimes'] = randTimes

   repeatSizeDif = margs.RepeatSizeDif; 
   if not repeatSizeDif==None:
      repeatSizeDifArray = repeatSizeDif.split('_'); usedRepeatSizeDifArray= []
      for irsda in range(len(repeatSizeDifArray)):
         if len(repeatSizeDifArray[irsda])>0: 
            usedRepeatSizeDifArray.append(int(repeatSizeDifArray[irsda]))
            errorStr += non_negative(usedRepeatSizeDifArray[-1], 'repeatSizeDif'+str(irsda));
      specifiedOptions['repeatSizeDif'] = usedRepeatSizeDifArray
      if len(usedRepeatSizeDifArray)>1:
         if usedRepeatSizeDifArray[0]>=usedRepeatSizeDifArray[1]:
            errorStr += ("\tError!!! the first repeat size difference (%s) must be smaller than the second (%s)" % (usedRepeatSizeDifArray[0], usedRepeatSizeDifArray[1]))+'\n'
   else:
      specifiedOptions['repeatSizeDif'] = []

   if not errorStr==originalError:
      print errorStr; #BAMinput|FASTQinput|Simulate
      parser.print_help();
      parser.parse_args(['Simulate', '--help']);
      sys.exit(140)

   simulation_file_id += ('_ins%.2f_del%.2f_sub%.2f_cov%d' % (insert_rate, del_rate, sub_rate, coverage))
   if not repeatSizeDif==None:
      simulation_file_id += ('_repeatSizeDif%s' % (repeatSizeDif))
   analysis_file_id += ('_times%d' % randTimes)
   if commonOptions['specifiedGeneName']==None:
      pass; #simulation_file_id = ("_def%s" % (commonOptions['repeatgene'])) + simulation_file_id
   else: simulation_file_id = ("_def%s" % (commonOptions['specifiedGeneName'])) + simulation_file_id

   unique_file_id =  simulation_file_id+analysis_file_id

   specifiedOptions['simulation_file_id'] = simulation_file_id
   specifiedOptions['analysis_file_id'] = analysis_file_id
   specifiedOptions['unique_file_id'] = unique_file_id

   #print commonOptions, specifiedOptions; #sys.exit(0);

   logfolder = 'logsim/'
   if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)
   filename = logfolder + 'TrinRepSim_' + commonOptions['repeatgene'] + unique_file_id + '.log'
   LOG_FILENAME = filename
   logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

   random.seed(7);

   if not commonOptions['specifiedGeneInfo'] == UserDefinedGenedefault:
        res = trinucleotideRepeatRealSimulation.getSimforKnownGeneWithPartialRev(commonOptions, specifiedOptions)
   else:
        if commonOptions['knownGeneName']=='all':
                res = trinucleotideRepeatRealSimulation.getSim(commonOptions, specifiedOptions);
        else:
                res = trinucleotideRepeatRealSimulation.getSimforKnownGene(commonOptions, specifiedOptions)

   for r1 in res:
        if not isUnsymAlign:
                logging.info('%s, %s<<<True<<<>>>Pred>>>%s, %s \t\t#wrongAlign %s' % (str(r1[0][0]), str(r1[0][1]), str(r1[1][0]), str(r1[1][1]), str(r1[2])));
        else: logging.info('%s, %s<<<True<<<>>>Pred>>>%s, %s bwa:%s %s' % (str(r1[0][0]), str(r1[0][1]), str(r1[1][0]), str(r1[1][1]), str(r1[2][0]), str(r1[2][1])));
        #print r1[0][0], r1[0][1], '<<<True<<<>>>Pred>>>', r1[1][0], r1[1][1], '\t\t#wrongAlign', r1[2];

   simresfolder = 'logsim/sim_res/'
   if not os.path.isdir(simresfolder):
        os.system('mkdir '+simresfolder)

   curfilename = simresfolder + commonOptions['repeatgene'] + unique_file_id + '.txt'
   curres = []
   for r1 in res:
        curstr = ('%s %s %s' % (str(r1[0][0]), str(r1[0][1]), str(r1[1][0])))
        if len(r1[1])>1 and (not string.strip(str(r1[1][1]))==''): curstr += (' %s' % str(r1[1][1]))
        else: curstr += (' %s' % str(r1[1][0]))
        if isUnsymAlign:
                curstr += (' %s %s' % (str(r1[2][0]), str(r1[2][1])))
        curres.append(curstr)
   findTrinucleotideRepeats.myWriteTxtFile(curres, curfilename);


def FASTQinput(margs):
   specifiedOptions = {}
   analysis_file_id = ""

   errorStr = originalError
   commonOptions, errorStr_com, analysis_file_id_com = getCommonOptions(margs)
   errorStr += errorStr_com

   isUnsymAlign = margs.UnsymAlign; specifiedOptions['isUnsymAlign'] = isUnsymAlign
   if isUnsymAlign: analysis_file_id += "_unsymalign"
   analysis_file_id += analysis_file_id_com

   specifiedOptions['fastafile'] = margs.fastq
   if specifiedOptions['fastafile']==None: errorStr += '\tNo fasta file provided\n';
   elif not os.path.isfile(specifiedOptions['fastafile']):  errorStr += '\tThe fasta file ('+specifiedOptions['fastafile']+') does not exit\n';

   if not errorStr==originalError:
      print errorStr; #BAMinput|FASTQinput|Simulate
      parser.print_help();
      parser.parse_args(['FASTQinput', '--help']);
      #ArgumentParser.print_help();
      sys.exit(140)

   unique_file_id =  commonOptions['specifiedGeneName']+analysis_file_id
   specifiedOptions['unique_file_id'] = unique_file_id
   specifiedOptions['analysis_file_id'] = analysis_file_id

   #print commonOptions, specifiedOptions; #sys.exit(0);

   logfolder = 'logfq/'
   if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)
   filename = logfolder + 'TrinRepFQ_' + commonOptions['repeatgene'] + unique_file_id + '.log'
   LOG_FILENAME = filename
   logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

   res = mySCA3.getSCA3forKnownGeneWithPartialRev(commonOptions, specifiedOptions)

   for r1 in res:
        print r1
        #print r1[:2], r1[2];
        pstr = '';
        if len(r1[0])>0: pstr += str(r1[0][0])
        if len(r1[0])>1: pstr += ','+str(r1[0][1])
        else: pstr += ','+str(r1[0][0])
        if isUnsymAlign and len(r1[1])>0:
           pstr += '('+str(r1[1][0])
           if len(r1[1])>1: pstr += str(r1[1][1])+')'
           else: pstr += str(r1[1][0])+')'
        logging.critical('\t Gene name='+commonOptions['repeatgene']+' repeat in test genome='+pstr)



def BAMinput(margs):
   specifiedOptions = {}
   analysis_file_id = ""

   errorStr = originalError
   commonOptions, errorStr_com, analysis_file_id_com = getCommonOptions(margs)
   errorStr += errorStr_com

   specifiedOptions["SepbamfileTemp"] = margs.SepbamfileTemp; 
   specifiedOptions["Onebamfile"] = margs.Onebamfile

   if not specifiedOptions["Onebamfile"]==None:
      specifiedOptions['bamfile'] = specifiedOptions["Onebamfile"]
      if not os.path.isfile(specifiedOptions['bamfile']):  errorStr += '\tThe bam file ('+specifiedOptions['bamfile']+') does not exit\n';
   else:
      specifiedOptions['bamfile'] = specifiedOptions["SepbamfileTemp"]
   if specifiedOptions['bamfile']==None: errorStr += '\tNo bam file provided\n';

   if not errorStr==originalError:
      print errorStr; #BAMinput|FASTQinput|Simulate
      parser.print_help();
      parser.parse_args(['BAMinput', '--help']);
      sys.exit(140)

   analysis_file_id += analysis_file_id_com
   if not specifiedOptions["Onebamfile"]==None:
      analysis_file_id += ("_for_%s" % (specifiedOptions['bamfile'].replace('/', '_')))
   else: 
      bamstr = specifiedOptions['bamfile'].replace('/', '_')
      bamstr = bamstr.replace('%s', '')
      analysis_file_id += ("_for_%s" % bamstr)

   unique_file_id = analysis_file_id
   specifiedOptions['analysis_file_id'] = analysis_file_id
   specifiedOptions['unique_file_id'] = unique_file_id

   #print (specifiedOptions['bamfile'] % '21')
   #print commonOptions, specifiedOptions; #sys.exit(0);

   logfolder = 'logbam/'
   if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)
   filename = logfolder + 'TrinRepBAM_' + commonOptions['repeatgene'] + unique_file_id + '.log'
   LOG_FILENAME = filename
   logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

   logging.info('Input BAM file is '+specifiedOptions['bamfile']);
   summary = []
   if margs.repeatgene=='all':
        summary = findTrinucleotideRepeats.getRepeat(commonOptions, specifiedOptions)
   else:
        summary.append(findTrinucleotideRepeats.getRepeatForKnownGene(commonOptions, specifiedOptions))

   for curgrep in summary:
       print curgrep[0], ('%.2f' % curgrep[1]), curgrep[2:]
       #print curgrep

   for curgrep in summary:
        pstr = 'NONE'
        if len(curgrep[2])>0: pstr = str(curgrep[2][0])
        if len(curgrep[2])>1: pstr += ','+ str(curgrep[2][1])
        else: 
            if len(curgrep[2])>0: pstr += ','+ str(curgrep[2][0])
            else: pstr += '0,0'
        logging.critical('\t Gene name='+curgrep[0]+'; ref_repeat='+('%.0f' % (curgrep[1]))+('(allreads=%d); ' % (curgrep[4]))+'repeat in test genome='+ pstr + ' #not_used_reads:'+str(curgrep[5]))
        logging.critical('\t#Format: Repeat_time=numberOfAlignedRead:'+curgrep[3])

   print ''
   gLoc = commonOptions['gLoc']; mprintstr = ''
   for curgrep in summary:
        pstr = 'NONE'
        if len(curgrep[2])>0: pstr = str(curgrep[2][0])
        if len(curgrep[2])>1: pstr += ','+ str(curgrep[2][1])
        else: 
            if len(curgrep[2])>0: pstr += ','+ str(curgrep[2][0])
            else: pstr += '0,0'
        logging.critical('\t Gene name='+curgrep[0]+' repeat in test genome='+ pstr)
   
        mprintstr += curgrep[0]+' '+pstr
        if gLoc.has_key(curgrep[0]) and len(gLoc[curgrep[0]][7])>0:
           mrange = trinucleotideRepeatRealSimulation.getRepeatRange2(gLoc[curgrep[0]][7])
           for mi in range(2):
              for ni in range(2):
                 if ni==0: mprintstr += ';' # ' ,\''
                 else: mprintstr += '-'
                 if len(mrange)>mi and len(mrange[mi])>ni: mprintstr += str(mrange[mi][ni]);
           
        mprintstr += '\n'
   print mprintstr


##############################################################################
#
#
#
##############################################################################

#parser = argparse.ArgumentParser(description='''Determine Trinucleotide repeat for (a) gene(s) of interests or for all genes.''', epilog="For example, \n \
#\tpython %(prog)s (BAMinput|FASTAQnput|Simulate); \n \
#\tpython %(prog)s (BAMinput|FASTQinput|Simulate) -repeatgene HTT; \n \
#\tpython %(prog)s (BAMinput|FASTQinput|Simulate) -repeatgene ATXN3 --correct 1 --updown 18 --extend 0; \n \
#\tpython %(prog)s (BAMinput|FASTQinput|Simulate) -repeatgene ATXN3 --correct 1 --updown 18 --extend 0 --coverage 100 --randTimes 100 -UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1; \n \
#OR \tpython %(prog)s -repeatgene all;\n \
#Final results was stored in log*/TrinRepSim*.log", formatter_class=RawTextHelpFormatter);


subparsers = parser.add_subparsers()
parent_parser = argparse.ArgumentParser(add_help=False)

########################################

com_group_for_align = parent_parser.add_argument_group('Common options for alignment')
com_group_for_align.add_argument("--hg", default='hg38', help="The reference genome is used. Currently, hg38 and hg19 are supported")
com_group_for_align.add_argument("--hgfile", default='', help="The file name of reference genome. It could be empty and the default file name is 'hg'.fa")
com_group_for_align.add_argument("--SeqDepth", default=0, help="The depth of the reads. For filtering in peak detection.")
com_group_for_align.add_argument("--RemInDel", type=int, default=1, help="Is unsymmetrical alignment used for error correction. 1: yes(Default), 0: no");
com_group_for_align.add_argument("--updown", type=int, default=18, help="Is upstream/downstream used for repeat length inference. non-0: yes(Default: 18), 0: no");
com_group_for_align.add_argument("--extend", type=int, default=3, help="Is upstream/downstream extended as repeat region. non-0: yes, 0: no(Default)");

####################################

com_group_for_gene = parent_parser.add_argument_group('Common options for gene information')
com_group_for_gene.add_argument("--repeatgene", default=None, help="A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. 'all': all known genes associated with trinucleotide repeat disorders will be analyzed;");
com_group_for_gene.add_argument("--UserDefinedGene", default=UserDefinedGenedefault, help="The gene information defined by users. If this option is given, the default gene information will be revised. Default: ///////");
com_group_for_gene.add_argument("--UserDefinedGeneName", default=None, help="The name for storing results. Default: Pcr1");

#####################################

parser_bam = subparsers.add_parser('BAMinput', parents=[parent_parser], help="Detect trinucleotide repeats from a BAM file", description="Detect trinucleotide repeats from a BAM file", epilog="For example, \n \
python %(prog)s --bamfile XXX.bam --repeatgene HTT; \n \
python %(prog)s --bamfile XXX.bam --repeatgene HTT --RemInDel 1 --updown 18 --extend 0; \n \
python %(prog)s --bamfile freeze4-all-merge.sort.bam --repeatgene HTT; \n \
Final results was stored in logbam/TrinRepBAM*.log \n \
Note!!!!!! \n \
BAM file should be produced with the same version of 'hg' if '--hg' or '--hgfile' is specified. \n \
", formatter_class=RawTextHelpFormatter)

bamgroup = parser_bam.add_mutually_exclusive_group(); #required=True)
bamgroup.add_argument("--Onebamfile", default=None, help="A BAM file storing all alignments")
bamgroup.add_argument("--SepbamfileTemp", default=None, help="A separated BAM file template storing all alignments; separated by chromosome ids. For example, '--SepbamfileTemp'=mybam_Chr%%s_sorted.bam where '%%s' is chromosome id (1....22,x,y)")

#parser_bam.add_argument("--bamfile", default=None, help="A BAM file storing all alignments");
parser_bam.set_defaults(func=BAMinput)

######################################

parser_fasta = subparsers.add_parser('FASTQinput', parents=[parent_parser], help="Detect trinucleotide repeats from a FASTQ file", description="Detect trinucleotide repeats from a FASTQ file", epilog="For example, \n \
python %(prog)s --fastq XXX.fq --repeatgene HTT --UnsymAlign 1; \n \
python %(prog)s --fastq XXX.fq --repeatgene HTT --RemInDel 1 --updown 18 --extend 0 --UnsymAlign 1; \n \
python %(prog)s --fastq XXX.fq --repeatgene HTT --RemInDel 1 --updown 18 --extend 0 --UnsymAlign 1 --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1; \n \
python %(prog)s --repeatgene atxn3 --UnsymAlign 1 --RemInDel 1 --updown 18 --extend 0 --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName sca3_pcr25_raw_test --fastq atxn3_data/rawdata/sam025.raw.fastq  \n \
Final results was stored in logfq/TrinRepFQ*.log \n \
", formatter_class=RawTextHelpFormatter)

parser_fasta.add_argument("--UnsymAlign", type=int, default=1, help="Whether unsymmetrical alignment (1) rather than bwa mem (0) would be used. Default: 1");
parser_fasta.add_argument("--fastq", help="The file name for fasta sequences");
parser_fasta.set_defaults(func=FASTQinput)

##############################################

parser_sim = subparsers.add_parser('Simulate', parents=[parent_parser], help="For simulation", description="Simulate for a given gene", epilog="For example, \n \
python %(prog)s --repeatgene HTT --UnsymAlign 1; \n \
python %(prog)s --repeatgene HTT --RemInDel 1 --updown 18 --extend 0 --UnsymAlign 1; \n \
python %(prog)s --repeatgene HTT --RemInDel 1 --updown 18 --extend 0 --UnsymAlign 1 --coverage 100 --randTimes 100 -UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1; \n \
python %(prog)s --repeatgene atxn3 --UnsymAlign 1 --RemInDel 1 --updown 18 --extend 0 --randTimes 100 --coverage 50 \n \
python %(prog)s --repeatgene atxn3 --UnsymAlign 1 --RemInDel 1 --updown 18 --extend 0 --randTimes 100 --coverage 50 --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName pcr1 \n \
Final results was stored in logsim/TrinRepSim*.log \n \
", formatter_class=RawTextHelpFormatter)

parser_sim.add_argument("--insert_rate", type=float, default=0.12, help="Insert error rate. Default: 0.12");
parser_sim.add_argument("--del_rate", type=float, default=0.02, help="Deletion error rate. Default: 0.02");
parser_sim.add_argument("--sub_rate", type=float, default=0.02, help="Substitution error rate. Default: 0.02");
parser_sim.add_argument("--coverage", type=int, default=300, help="The number of reads produced after mutations. Default: 300");
parser_sim.add_argument("--randTimes", type=int, default=100, help="The number of simulation times. Default: 100");
parser_sim.add_argument("--UnsymAlign", type=int, default=1, help="Whether unsymmetrical alignment (1) rather than bwa mem (0) would be used. Default: 1");
parser_sim.add_argument("--RepeatSizeDif", default=None, help="The absolute difference of repeat size in two alleles: Default: None; could be 1, 2, or 3_4 and so on. Must be positive numbers.");
parser_sim.set_defaults(func=simulation)


#parser.print_help();

if len(sys.argv)<2:
   parser.print_help();
else:
   args = parser.parse_args()
   args.func(args);



