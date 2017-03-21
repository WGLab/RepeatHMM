
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
from scripts import mySCA3
from scripts.myheader import *
from scripts import printHMMmatrix



parser = argparse.ArgumentParser(description="Determine microsatellite repeat for (a) gene(s) of interests or for all genes.", epilog="For example, \n \
\tpython %(prog)s BAMinput: with a BAM file as input\n \
\tpython %(prog)s FASTQinput: with a FASTQ file as input \n \
\tpython %(prog)s Simulate: for simulation", formatter_class=RawTextHelpFormatter);


originalError = '!!!Error: !!!!!! \n'


def non_negative(i, mstr):
   if i<0: return (("\tError %d could not be negative(%d)\n" % (mstr, i)))
   else: return ''

def check_TRFOptions(trfo):
	trfosp = trfo.split('_')
	errorMsg = ['']
	if len(trfosp)==7 or len(trfosp)==6: 
		mstrs = ['Match(1st)', 'Mismatch(2rd)', 'Delta(3rd)', 'PM(4th)', 'PI(5th)', 'Minscore(6th)', 'MaxPeriod(7th)']
		for imst in range(len(mstrs)):
			if len(trfosp)==6 and imst==6: continue;
			retmsg = non_negative(int(trfosp[imst]), mstrs[imst])
			if not retmsg=='': errorMsg.append.append(retmsg)
	else:
		errorMsg.append("The number of TRFOptions is not 7:")
		errorMsg.append(str(len(trfosp)))
		errorMsg.append('\n')

	return ''.join(errorMsg)

def checkM(mm, colnum=None):
   yourrows = mm.split(';'); expnum = 100;
   yourm = []; rowsize = {}; sumsize = {}
   for i in range(len(yourrows)):
      yourcr = []
      rv = yourrows[i].split(',');
      cursum = 0;
      for curv in rv:
         yourcr.append(float(curv));
         cursum += yourcr[-1]*expnum
         if yourcr[-1]<1e-9: yourcr[-1] = 1e-9
      yourm.append(yourcr);
      if not rowsize.has_key(len(yourcr)):
         rowsize[len(yourcr)] = []
      rowsize[len(yourcr)].append(i)
      cursum = round(int(cursum+0.5)/float(expnum), 3)
      if not sumsize.has_key(cursum):
         sumsize[cursum] = []
      sumsize[cursum].append(i)

   errormsg = []; wrongrowsize=False;
   if len(rowsize)>1: 
      errormsg.append('Error: size for each row is not the same:'); wrongrowsize=True;
   elif len(rowsize)<1:
      errormsg.append('Error: size for each row is less than 1:'); wrongrowsize=True;
   else:
      if (not colnum==None):
         if rowsize.keys()[0]==colnum: pass
         else: 
            errormsg.append('Error: the number of columns is not 4'); wrongrowsize=True;
      else:
         if rowsize.keys()[0]==len(yourm): pass
         else:
            errormsg.append('Error: the number of columns (%d) is not equal to the number of rows (%d)' % (rowsize.keys()[0], len(yourm))); wrongrowsize=True;
   if wrongrowsize: errormsg.append(str(rowsize))
   
   if (not len(sumsize)==1) or (len(sumsize)==1 and not (expnum-0.1<sumsize.keys()[0]<expnum+0.1)):
      errormsg.append('Error: the sum of each row is not correct:'); errormsg.append(str(sumsize))

   return yourm, '\n'.join(errormsg)

def setInsDelSub(comOpt):
   #["Pacbio", "Nanopore", "Illumina", None]
   if comOpt['SeqTech']=="Pacbio":
      comOpt['hmm_insert_rate'], comOpt['hmm_del_rate'], comOpt['hmm_sub_rate'] = 0.110, 0.020, 0.020
   elif comOpt['SeqTech']=="Nanopore":
      comOpt['hmm_insert_rate'], comOpt['hmm_del_rate'], comOpt['hmm_sub_rate'] = 0.100, 0.050, 0.050
   elif comOpt['SeqTech']=="Illumina":
      comOpt['hmm_insert_rate'], comOpt['hmm_del_rate'], comOpt['hmm_sub_rate'] = 0.002, 0.002, 0.002

def getCommonOptions(margs):
   errorStr = ""
  
   isRemInDel, isupdown, isExtend = margs.RemInDel, margs.updown, margs.extend
   errorStr += non_negative(isRemInDel, 'isRemInDel');
   errorStr += non_negative(isupdown, 'isupdown');
   errorStr += non_negative(isExtend, 'isExtend');

   MinSup = margs.MinSup;
   errorStr += non_negative(MinSup, 'MinSup');

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

   SplitAndReAlign = margs.SplitAndReAlign
   TRFOptions = margs.TRFOptions
   if not SplitAndReAlign==0:
      errorStr += check_TRFOptions(TRFOptions)
      analysis_file_id_common += '_SplitAndReAlign'
      analysis_file_id_common += '_'+TRFOptions

   minTailSize = margs.minTailSize
   if minTailSize<10: minTailSize = 10;
   minRepBWTSize = margs.minRepBWTSize
   if minRepBWTSize<10: minRepBWTSize=10;
   RepeatTime = margs.RepeatTime
   if RepeatTime<2: RepeatTime = 2
 
   BWAMEMOptions = margs.BWAMEMOptions.replace('_', ' -')
   if not BWAMEMOptions[0]=='-':  BWAMEMOptions = ' -' + BWAMEMOptions
   if not BWAMEMOptions[-1]==' ': BWAMEMOptions = BWAMEMOptions + ' '
   MaxRep = margs.MaxRep
   if MaxRep<14: MaxRep = 14;
   CompRep = printHMMmatrix.getCompRep(margs.CompRep)

   commonOptions = {}
   commonOptions['MinSup'] = MinSup
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

   commonOptions['SplitAndReAlign'] = SplitAndReAlign
   commonOptions['TRFOptions'] = TRFOptions
   commonOptions['minTailSize'] = minTailSize
   commonOptions['minRepBWTSize'] = minRepBWTSize
   commonOptions['RepeatTime'] = RepeatTime

   commonOptions['BWAMEMOptions'] = BWAMEMOptions
   commonOptions['MaxRep'] = MaxRep
   commonOptions['CompRep'] = CompRep

   analysis_file_id_common += '_'+commonOptions['hg']+'_comp'

   if not commonOptions['specifiedGeneName']==None:
      analysis_file_id_common += '_def'+commonOptions['specifiedGeneName']

   hmm_insert_rate = margs.hmm_insert_rate
   hmm_del_rate = margs.hmm_del_rate
   hmm_sub_rate = margs.hmm_sub_rate
   SeqTech = margs.SeqTech
   transitionm = margs.transitionm
   emissionm = margs.emissionm
   commonOptions['hmm_insert_rate'] = hmm_insert_rate
   commonOptions['hmm_del_rate'] = hmm_del_rate
   commonOptions['hmm_sub_rate'] = hmm_sub_rate
   commonOptions['SeqTech'] = SeqTech
   if not SeqTech==None: 
      analysis_file_id_common += '_'+SeqTech
      setInsDelSub(commonOptions)

   analysis_file_id_common += ('_I%.3f' % hmm_insert_rate)
   analysis_file_id_common += ('_D%.3f' % hmm_del_rate)
   analysis_file_id_common += ('_S%.3f' % hmm_sub_rate)

   if not transitionm==None: 
      commonOptions['transitionm'], cerrstr = checkM(commonOptions['transitionm'], None)
      if not cerrstr=="": errorStr += ''.join(['\t', cerrstr, '\n'])
   else: commonOptions['transitionm'] = transitionm
   if not emissionm==None: 
      commonOptions['emissionm'], cerrstr = checkM(commonOptions['emissionm'], 4)
      if not cerrstr=="": errorStr += ''.join(['\t', cerrstr, '\n'])
   else: commonOptions['emissionm'] = emissionm

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

   specifiedOptions['simulation_file_id'] = commonOptions['repeatgene'] + simulation_file_id
   specifiedOptions['analysis_file_id'] = commonOptions['repeatgene'] + analysis_file_id
   specifiedOptions['unique_file_id'] = commonOptions['repeatgene'] + unique_file_id

   logfolder = 'logsim/'
   if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)
   filename = logfolder + 'TrinRepSim_' + unique_file_id + '.log'
   LOG_FILENAME = filename
   logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

   random.seed(7);

   if not commonOptions['specifiedGeneInfo'] == UserDefinedGenedefault:
      summary= trinucleotideRepeatRealSimulation.getSimforKnownGeneWithPartialRev(commonOptions, specifiedOptions)
   else:
      summary = trinucleotideRepeatRealSimulation.getSimforKnownGene(commonOptions, specifiedOptions)

   printRepInfo(summary)

   simresfolder = 'logsim/sim_res/'
   if not os.path.isdir(simresfolder):
      os.system('mkdir '+simresfolder)

   curfilename = simresfolder + commonOptions['repeatgene'] + unique_file_id + '_comp.txt'
   allres = []
   summarykeys = summary.keys(); summarykeys.sort()
   for sumk in summarykeys:
      if sumk==summarykeys[0]:
         titstr = ['#No']
      methtypekeys = summary[sumk][0].keys(); methtypekeys.sort();
      curres = [str(sumk)]
      for methk in methtypekeys:
         if sumk==summarykeys[0]:
            titstr.append(methk)
         res1 = summary[sumk][0][methk]
         curres.append('%s %s' % (str(res1[0]), str(res1[1])))
      if sumk==summarykeys[0]:
         allres.append(' '.join(titstr))
      allres.append(' '.join(curres))
		
   findTrinucleotideRepeats.myWriteTxtFile(allres, curfilename);


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
      sys.exit(140)

   unique_file_id =  commonOptions['specifiedGeneName']+analysis_file_id
   specifiedOptions['unique_file_id'] = unique_file_id
   specifiedOptions['analysis_file_id'] = analysis_file_id

   logfolder = 'logfq/'
   if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)
   filename = logfolder + 'TrinRepFQ_' + commonOptions['repeatgene'] + unique_file_id + '.log'
   LOG_FILENAME = filename
   logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

   summary = {}
   summary[commonOptions['repeatgene']] = mySCA3.getSCA3forKnownGeneWithPartialRev(commonOptions, specifiedOptions)
   print '\nfor output'; logging.info('\nfor output')
   printRepInfo(summary)


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
      pass #analysis_file_id += ("_for_%s" % (specifiedOptions['bamfile'].replace('/', '_')))
   else: 
      bamstr = specifiedOptions['bamfile'].replace('/', '_')
      bamstr = bamstr.replace('%s', '')
      #analysis_file_id += ("_for_%s" % bamstr)

   unique_file_id = analysis_file_id
   specifiedOptions['analysis_file_id'] = analysis_file_id
   specifiedOptions['unique_file_id'] = unique_file_id

   logfolder = 'logbam/'
   if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)
   filename = logfolder + 'TrinRepBAM_' + commonOptions['repeatgene'] + unique_file_id + '.log'
   LOG_FILENAME = filename
   logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

   logging.info('Input BAM file is '+specifiedOptions['bamfile']);

   if margs.repeatgene=='all':
      summary = findTrinucleotideRepeats.getRepeat(commonOptions, specifiedOptions)
   else:
      summary = {}
      summary[margs.repeatgene] = findTrinucleotideRepeats.getRepeatForKnownGene(commonOptions, specifiedOptions)

   printRepInfo(summary)

def printRepInfo(summary):
	print ''; logging.info('')
	sumkeys = summary.keys(); sumkeys.sort()
	for sumk in sumkeys:
		methkeys = summary[sumk][1].keys(); methkeys.sort()
		for mk in methkeys:
			print summary[sumk][1][mk]
			logging.info(summary[sumk][1][mk]+'')

	print ''; logging.info('')
	print ''; logging.info('')
	print ''; logging.info('')

	for sumk in sumkeys:
		if sumk==sumkeys[0]:
			titstr = [' ']

		prstr = [str(sumk)]
		methkeys = summary[sumk][0].keys(); methkeys.sort()
		for mk in methkeys:
			rep1 = summary[sumk][0][mk]
			prstr.append(' %d %d;' % (rep1[0], rep1[1])) 
			if sumk==sumkeys[0]:
				titstr.append(mk)
		if sumk==sumkeys[0]:
			print '\t'.join(titstr)
			logging.info(('\t'.join(titstr))+'')
		print '\t'.join(prstr)
		logging.info(('\t'.join(prstr))+'')
	print ''; logging.info('')	


##############################################################################
#
#
#
##############################################################################



subparsers = parser.add_subparsers()
parent_parser = argparse.ArgumentParser(add_help=False)

########################################

com_group_for_align = parent_parser.add_argument_group('Common options for alignment')
com_group_for_align.add_argument("--hg", default='hg38', help="The reference genome is used. Currently, hg38 and hg19 are supported")
com_group_for_align.add_argument("--hgfile", default='', help="The file name of reference genome. It could be empty and the default file name is 'hg'.fa")
com_group_for_align.add_argument("--RemInDel", type=int, default=1, help="Is unsymmetrical alignment used for error correction. 1: yes(Default), 0: no");
com_group_for_align.add_argument("--updown", type=int, default=18, help="Is upstream/downstream used for repeat length inference. non-0: yes(Default: 18), 0: no");
com_group_for_align.add_argument("--extend", type=int, default=3, help="Is upstream/downstream extended as repeat region. non-0: yes, 0: no(Default)");

##################################--BWAMEMOptions --MaxRep --CompRep
com_group_for_others = parent_parser.add_argument_group('Common options')
com_group_for_others.add_argument("--MinSup", type=int, default=5, help="The minimum reads associated with peaks.")
com_group_for_others.add_argument("--MaxRep", type=int, default=4000, help="The maximum repeat size. The smallest MaxRep should not be less than 14.")
com_group_for_others.add_argument("--CompRep", default='0', help="Whether the repeat pattern is simple ('0') or complex (nonzeo-'0': AlTlT50/C50lClT/C). For complex patterns, all patterns are required to have the same length, each position is seperated by `l` and the nucleotides at the same position is seperated by '/' where a nucleotide can by followed by a number specify relative frequency (cannot be float).")

####################################

com_group_for_gene = parent_parser.add_argument_group('Common options for gene information')
com_group_for_gene.add_argument("--repeatgene", default=None, help="A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. 'all': all known genes associated with trinucleotide repeat disorders will be analyzed;");
com_group_for_gene.add_argument("--UserDefinedGene", default=UserDefinedGenedefault, help="The gene information defined by users. If this option is given, the default gene information will be revised. Default: ///////");
com_group_for_gene.add_argument("--UserDefinedGeneName", default=None, help="The name for storing results. Default: Pcr1");

com_group_for_splitalign = parent_parser.add_argument_group('Common options for re-alignment after splitting long reads using repeat regions')
com_group_for_splitalign.add_argument("--SplitAndReAlign", type=int, choices=[0, 1, 2], default=0, help="Split long reads using repeat region in long reads and re-align the non-repeat regions using BWA MEM. Default=0: not use SplitAndReAlign; 1: use SplitAndReAlign only; 2: combine 0 and 1")
com_group_for_splitalign.add_argument("--TRFOptions",  default="2_7_4_80_10_100", help="The options used for detecting repeat region in a read using Tandem Repeat Finder. The options are merging using _. Default='2_7_4_80_10_100_500' or '2_7_4_80_10_100'. The last parameter will be twice of the length of repeat if not given.")
com_group_for_splitalign.add_argument("--minTailSize", type=int, default=70, help="After the split using repeat regions, discard the leftmost/rightmost non-repeat sub-sequences if they have less than --minTailSize bps. Must not be less than 10. Default=70")
com_group_for_splitalign.add_argument("--minRepBWTSize", type=int, default=70, help="After the split using repeat regions, merge any two non-repeat sub-sequences if they distance is less than --minRepBWTSize bps. Must not be less than 10. Default=70")
com_group_for_splitalign.add_argument("--RepeatTime", type=int, default=5, help="The minimum repeat time for a microsatellite detected by Tandem Repeat Finder. Must not be less than 2. Default=5")
com_group_for_splitalign.add_argument("--BWAMEMOptions",  default="k8_W8_r7", help="The options used for BWA MEM to align sub-reads after splitting. The options are merging using _. Default='k8_W8_r7'. For example: 'k8_W8_r7'")

#fafqfile\|fafqtype
#hmm_insert_rate hmm_del_rate hmm_sub_rate SeqTech transitionm emissionm
com_group_for_hmm = parent_parser.add_argument_group('Common options for setting HMM matrices')
com_group_for_hmm.add_argument("--hmm_insert_rate", type=float, default=0.12, help="Insert error rate in long reads. Default: 0.12");
com_group_for_hmm.add_argument("--hmm_del_rate", type=float, default=0.02, help="Deletion error rate in long reads. Default: 0.02");
com_group_for_hmm.add_argument("--hmm_sub_rate", type=float, default=0.02, help="Substitution error rate in long reads. Default: 0.02");
com_group_for_hmm.add_argument("--SeqTech", choices=["Pacbio", "Nanopore", "Illumina", None], default=None, help="The sequencing techniques, Pacbio or Nanopore or Illumina, used to generate reads data. Default: None. Setting this option will override the setting for --hmm_insert_rate, --hmm_del_rate and --hmm_sub_rate");
com_group_for_hmm.add_argument("--transitionm", default=None, help="User-specified transition matrix for HMM. The number of rows and columns must be the same as the 3*L+1 where L is the size of repeat unit. Probabilities in a row is separated by ',' and their sum must be 1. Probabilities of different rows are separated by ';'. Please pay more attention when providing this parameter. For CG repeat, the example of this matrix is '0.96,0.02,0,0,0,0.02,0;0,0.001,0.869,0.11,0,0,0.02;0.02,0.849,0.001,0,0.11,0.02,0;0,0.001,0.869,0.11,0,0,0.02;0.02,0.849,0.001,0,0.11,0.02,0;0.02,0.849,0.001,0,0.11,0.02,0;0.02,0.001,0.849,0.11,0,0,0.02'. Setting this option will override the setting of --hmm_insert_rate, --hmm_del_rate, --hmm_sub_rate and --SeqTech for transition matrix.");
com_group_for_hmm.add_argument("--emissionm", default=None, help="User-specified emission matrix for HMM. The number of rows must be the same as the 3*L+1 where L is the size of repeat unit and the number of columns must be 4. Probabilities in a row is separated by ',' and their sum must be 1. Probabilities of different rows are separated by ';'. Please pay more attention when providing this parameter. For CG repeat, the example of this matrix is '0.25,0.25,0.25,0.25;0.005,0.985,0.005,0.005;0.005,0.005,0.985,0.005;0.25,0.25,0.25,0.25;0.25,0.25,0.25,0.25;0.005,0.005,0.985,0.005;0.005,0.985,0.005,0.005'. Setting this option will override the setting of --hmm_insert_rate, --hmm_del_rate, --hmm_sub_rate and --SeqTech for emission matrix.");

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


if len(sys.argv)<2:
   parser.print_help();
else:
   args = parser.parse_args()
   args.func(args);





