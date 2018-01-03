
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

from scripts import myBAMhandler #
from scripts import myFASTQhandler
from scripts.myheader import *
from scripts import printHMMmatrix
from scripts import myPredefinedPatternReader
from scripts import myScanWholeGenome
from scripts import myCommonFun


parser = argparse.ArgumentParser(description="Determine microsatellite repeat of interests or for all microsatellites.", epilog="For example, \n \
\tpython %(prog)s BAMinput: with a BAM file as input\n \
\tpython %(prog)s FASTQinput: with a FASTQ file as input \n \
", formatter_class=RawTextHelpFormatter);
#\tpython %(prog)s Scan: for Scaning whole genome", formatter_class=RawTextHelpFormatter);


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

def printOptions(op):
	opkeys = op.keys(); opkeys.sort();
	for opk in opkeys:
		if (opk not in ['gLoc']): #(not op[opk]==None) and (opk not in ['gLoc']):
			print ''.join([('%20s' % opk), '\t(', str(op[opk]), ');'])
	print ''

def getCommonOptions(margs, cominfo=None):
   random.seed(7);
   commonOptions = {}
   errorStr = ""
 
   if margs.outlog<0:
      commonOptions['outlog'] = M_WARNING
   else: commonOptions['outlog'] = margs.outlog
 
   isGapCorrection, repeatFlankLength = margs.GapCorrection, margs.FlankLength
   errorStr += non_negative(isGapCorrection, 'isGapCorrection');
   errorStr += non_negative(repeatFlankLength, 'repeatFlankLength');

   MinSup = margs.MinSup;
   errorStr += non_negative(MinSup, 'MinSup');

   MatchInfo = margs.MatchInfo
   if MatchInfo==None:
      #commonOptions['MatchInfo'] = [3, -2, -2, -15, -1]
      commonOptions['MatchInfo'] = [2, -1, -2, -13, -1]
   else:
      mi_sp = MatchInfo.split(';')
      for i in range(len(mi_sp)):
         mi_sp[i] = int(mi_sp[i])
         if i==0:
            errorStr += non_negative(mi_sp[i], 'MatchInfo['+str(i+1)+']');
         else:
            errorStr += non_negative(-mi_sp[i], 'MatchInfo['+str(i+1)+']');
      commonOptions['MatchInfo'] = mi_sp

   analysis_file_id_common = ''
   analysis_file_id_common += '_GapCorrection' + str(isGapCorrection)
   analysis_file_id_common += '_FlankLength'+str(repeatFlankLength)

   commonOptions['hg'] = margs.hg
   if margs.hgfile==None or margs.hgfile=='':
      #commonOptions['hgfile'] = commonOptions['hg']+'.fa'
      #if not os.path.isfile(hg_reference_and_index+'/'+margs.hg+'.fa'):
      #   errorStr += '\tNo reference genome under the folder '+hg_reference_and_index+' with the filename '+margs.hg+'.fa\n'
      #else: 
      commonOptions['hgfile'] = hg_reference_and_index+'/'+margs.hg+'.fa'
   else:
      commonOptions['hgfile'] = margs.hgfile
   if not os.path.isfile(commonOptions['hgfile']):
      errorStr += '\tA reference file is needed (specified by "--hgfile" or under the folder '+hg_reference_and_index+' with the filename '+margs.hg+'.fa\n'

   #if not os.path.isfile(hg_reference_and_index+'/'+margs.hg+'.fa'):
   #   errorStr += '\tNo reference genome under the folder '+hg_reference_and_index+' with the filename '+margs.hg+'.fa\n'
   
   UserDefinedUniqID = margs.UserDefinedUniqID
   RepeatName = margs.repeatName
   specifiedRepeatInfo = margs.UserDefinedRepeat

   if RepeatName==None and (cominfo==None or ((not cominfo==None) and (cominfo.has_key('scan') and not cominfo['scan']==1))):
      errorStr += '\tNone gene/repeat information is given. \n'
   if (not specifiedRepeatInfo==UserDefinedRepeatdefault) and RepeatName=='all':
      errorStr += '\t"repeatName" cannot be "all" when specifying "UserDefinedRepeat"\n'

   SplitAndReAlign = margs.SplitAndReAlign
   TRFOptions = margs.TRFOptions
   if not SplitAndReAlign==0:
      errorStr += check_TRFOptions(TRFOptions)
      analysis_file_id_common += '_SplitAndReAlign'+str(SplitAndReAlign)
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

   commonOptions['MinSup'] = MinSup
   commonOptions['isGapCorrection'] = isGapCorrection
   commonOptions['repeatFlankLength'] = repeatFlankLength

   commonOptions['UserDefinedUniqID'] = UserDefinedUniqID
   commonOptions['repeatName'] = RepeatName
   commonOptions['specifiedRepeatInfo'] = specifiedRepeatInfo
   #commonOptions['hg'] = margs.hg
   #if margs.hgfile==None or margs.hgfile=='':
   #   commonOptions['hgfile'] = commonOptions['hg']+'.fa'
   #else:
   #   commonOptions['hgfile'] = margs.hgfile
   #if not os.path.isfile(commonOptions['hgfile']):
   #   errorStr += '\tA reference file is needed (specified by "--hgfile").\n'

   commonOptions['SplitAndReAlign'] = SplitAndReAlign
   commonOptions['TRFOptions'] = TRFOptions
   commonOptions['minTailSize'] = minTailSize
   commonOptions['minRepBWTSize'] = minRepBWTSize
   commonOptions['RepeatTime'] = RepeatTime

   commonOptions['BWAMEMOptions'] = BWAMEMOptions
   commonOptions['MaxRep'] = MaxRep
   commonOptions['CompRep'] = CompRep
   commonOptions['Tolerate_mismatch'] = margs.Tolerate

   analysis_file_id_common += '_'+commonOptions['hg']+'_comp'

   if not commonOptions['UserDefinedUniqID']==None:
      commonOptions['UserDefinedUniqID'] = commonOptions['UserDefinedUniqID'].replace(':', '_')
      commonOptions['UserDefinedUniqID'] = commonOptions['UserDefinedUniqID'].replace('/', '_')
      commonOptions['UserDefinedUniqID'] = commonOptions['UserDefinedUniqID'].replace('\\', '_')
      commonOptions['UserDefinedUniqID'] = commonOptions['UserDefinedUniqID'].replace('-', '_')
      analysis_file_id_common += '_'+commonOptions['UserDefinedUniqID']

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

   commonOptions['stsBasedFolder'] = stsBasedFolder
   moptions = {}
   Patternfile = margs.Patternfile
   if not Patternfile==None:
      commonOptions['Patternfile'] = Patternfile.split(';')
      for fi in commonOptions['Patternfile']:
         if fi[-3:]=='.pa':
            moptions['pafile'] = fi
   else: commonOptions['Patternfile'] = None

   moptions['stsBasedFolder'] = stsBasedFolder
   moptions['hg'] = margs.hg
   commonOptions['gLoc'] = myPredefinedPatternReader.getPredefinedMicrosatellites(moptions)

   if not margs.outFolder[-1]=='/': 
      commonOptions['align'] = margs.outFolder + '/'
   else: commonOptions['align'] = margs.outFolder
   if not os.path.isdir(margs.outFolder):
      os.system('mkdir '+margs.outFolder)

   return commonOptions, errorStr, analysis_file_id_common

def scan(margs):
   print 'This function is disenabled now.'
   print 'We will publish this function later.'
   sys.exit(140)

   cominfo = {}
   cominfo['scan'] = 1;
   #print margs;
   specifiedOptions = {}
   analysis_file_id = ""

   errorStr = originalError
   commonOptions, errorStr_com, analysis_file_id_com = getCommonOptions(margs, cominfo)
   errorStr += errorStr_com
  
   analysis_file_id += analysis_file_id_com

   if not errorStr==originalError:
      print errorStr; #BAMinput|FASTQinput|Scan
      parser.print_help();
      parser.parse_args(['Scan', '--help']);
      sys.exit(140)

   scan_region = margs.region;
   specifiedOptions['scan_region'] = scan_region
   if scan_region==None:
      specifiedOptions['analysis_file_id'] = 'gmm' + analysis_file_id
   else:
      scan_region_rp = scan_region.replace(':','_')
      scan_region_rp = scan_region_rp.replace('-','_')
      if analysis_file_id.find(scan_region_rp)==-1:		
         specifiedOptions['analysis_file_id'] = 'gmm_'+ scan_region_rp + analysis_file_id
      else: specifiedOptions['analysis_file_id'] = 'gmm'+ analysis_file_id

   specifiedOptions['unique_file_id'] = specifiedOptions['analysis_file_id']

   specifiedOptions["SepbamfileTemp"] = margs.SepbamfileTemp;
   specifiedOptions["Onebamfile"] = margs.Onebamfile

   if not specifiedOptions["Onebamfile"]==None:
      specifiedOptions['bamfile'] = specifiedOptions["Onebamfile"]
      if not os.path.isfile(specifiedOptions['bamfile']):  errorStr += '\tThe bam file ('+specifiedOptions['bamfile']+') does not exit\n';
   else:
      specifiedOptions['bamfile'] = specifiedOptions["SepbamfileTemp"]
   if specifiedOptions['bamfile']==None: errorStr += '\tNo bam file provided\n';

   specifiedOptions['align'] = commonOptions['align']
   #specifiedOptions['align'] = margs.outFolder
   #if not os.path.isdir(margs.outFolder):
   #   os.system('mkdir '+margs.outFolder)
   commonOptions['outlog'] = M_WARNING
   logfolder = 'logscan/'
   if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)
   filename = logfolder + 'RepScan_' + specifiedOptions['analysis_file_id'] + '.log'
   LOG_FILENAME = filename
   logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

   print 'The following options are used (included default):'
   printOptions(commonOptions)
   printOptions(specifiedOptions)

   specifiedOptions['thread'] = margs.thread
   specifiedOptions['continue'] = margs.conted
   specifiedOptions['scanresfolder'] = 'logscan/scan_res/'
   if margs.thread==1:
      mres, mdetail = myScanWholeGenome.scan(commonOptions, specifiedOptions)
      myCommonFun.myWriteScanResults(specifiedOptions, mres, mdetail, procss_info='')
   else:
      mres, mdetail = myScanWholeGenome.scan_multiprocess(commonOptions, specifiedOptions)
      myCommonFun.myWriteScanResults(specifiedOptions, mres, mdetail, procss_info='_all')

def FASTQinput(margs):
   specifiedOptions = {}
   analysis_file_id = ""

   errorStr = originalError
   commonOptions, errorStr_com, analysis_file_id_com = getCommonOptions(margs)
   errorStr += errorStr_com

   analysis_file_id += analysis_file_id_com

   specifiedOptions['fastafile'] = margs.fastq
   if specifiedOptions['fastafile']==None: errorStr += '\tNo fasta file provided\n';
   elif not os.path.isfile(specifiedOptions['fastafile']):  errorStr += '\tThe fasta file ('+specifiedOptions['fastafile']+') does not exit\n';

   if not errorStr==originalError:
      print errorStr; #BAMinput|FASTQinput|Scan
      parser.print_help();
      parser.parse_args(['FASTQinput', '--help']);
      sys.exit(140)

   if not commonOptions['UserDefinedUniqID']==None:
      unique_file_id =  commonOptions['UserDefinedUniqID']+analysis_file_id
   else: unique_file_id = analysis_file_id
   unique_file_id = '.gmm' + unique_file_id
   specifiedOptions['unique_file_id'] = unique_file_id
   specifiedOptions['analysis_file_id'] = analysis_file_id

   logfolder = 'logfq/'
   if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)
   filename = logfolder + 'RepFQ_' + commonOptions['repeatName'] + unique_file_id + '.log'
   LOG_FILENAME = filename
   logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

   specifiedOptions['align'] = commonOptions['align']

   print 'The following options are used (included default):'
   printOptions(commonOptions)
   printOptions(specifiedOptions)
   summary = {}
   summary[commonOptions['repeatName']] = myFASTQhandler.getSCA3forKnownGeneWithPartialRev(commonOptions, specifiedOptions)
   print '\nfor output'; logging.info('\nfor output')
   printRepInfo(summary)

   print 'The result is in', filename

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
      print errorStr; #BAMinput|FASTQinput|Scan
      parser.print_help();
      parser.parse_args(['BAMinput', '--help']);
      sys.exit(140)

   analysis_file_id += analysis_file_id_com

   unique_file_id = '.gmm' + analysis_file_id
   specifiedOptions['analysis_file_id'] = analysis_file_id
   specifiedOptions['unique_file_id'] = unique_file_id

   logfolder = 'logbam/'
   if not os.path.isdir(logfolder):
        os.system('mkdir '+logfolder)
   filename = logfolder + 'RepBAM_' + commonOptions['repeatName'] + unique_file_id + '.log'
   LOG_FILENAME = filename
   logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

   logging.info('Input BAM file is '+specifiedOptions['bamfile']);

   specifiedOptions['align'] = commonOptions['align']

   print 'The following options are used (included default):'
   printOptions(commonOptions)
   printOptions(specifiedOptions)

   if margs.repeatName=='all':
      summary = myBAMhandler.getRepeat(commonOptions, specifiedOptions)
   else:
      summary = {}
      summary[margs.repeatName] = myBAMhandler.getRepeatForKnownGene(commonOptions, specifiedOptions)

   printRepInfo(summary)

   print 'The result is in', filename

def printRepInfo(summary):
	print ''; logging.info('')
	sumkeys = summary.keys(); sumkeys.sort()
	for sumk in sumkeys:
		methkeys = summary[sumk][1].keys(); methkeys.sort()
		for mk in methkeys:
			print summary[sumk][1][mk]
			logging.info(summary[sumk][1][mk]+'')

	print ''; logging.info('')
	#print ''; logging.info('')
	#print ''; logging.info('')

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
# --hg --hgfile --GapCorrection --FlankLength --MatchInfo
com_group_for_align = parent_parser.add_argument_group('Common options for alignment')
com_group_for_align.add_argument("--hg", default='hg38', help="The reference genome is used. Currently, some repeat information in hg38 and hg19 are pre-defined. The default folder is ./mhgversion/")
com_group_for_align.add_argument("--hgfile", default='', help="The file name of reference genome. It could be empty and the default file name is 'hg'.fa")
com_group_for_align.add_argument("--GapCorrection", type=int, default=1, help="Is unsymmetrical alignment used for error correction of detected repeat region in reads. 1: yes(Default), 0: no");
com_group_for_align.add_argument("--FlankLength", type=int, default=30, help="Is flanking sequence used for repeat region detection. non-0: yes(Default: 30), 0: no");
com_group_for_align.add_argument("--MatchInfo", default="3;-2;-2;-15;-1", help="The match strategy for gap correction.");

################################## --Tolerate --MinSup --MaxRep --CompRep
com_group_for_others = parent_parser.add_argument_group('Common options')
com_group_for_others.add_argument("--outlog", default=M_WARNING, help="The level for output of running information: 0:DEBUG; 1:INFO; 2:WARNING; 3:ERROR; 4:FATAL")
com_group_for_others.add_argument("--Tolerate",  default=None, help="Tolerate mismatch, e.g., CTG:TTG:0:2 or CTG:CTT:-2:0.")
com_group_for_others.add_argument("--MinSup", type=int, default=5, help="The minimum reads associated with peaks.")
com_group_for_others.add_argument("--MaxRep", type=int, default=4000, help="The maximum repeat size. The smallest MaxRep should not be less than 14.")
com_group_for_others.add_argument("--CompRep", default='0', help="Whether the repeat pattern is simple ('0') or complex (nonzeo-'0': AlTlT50/C50lClT/C). For complex patterns, all patterns are required to have the same length, each position is seperated by `l` and the nucleotides at the same position is seperated by '/' where a nucleotide can by followed by a number specify relative frequency (cannot be float).")

com_group_for_others.add_argument("--outFolder", default='align/', help="The default folder for temporary output. Default: align/");

#################################### --repeatName --UserDefinedUniqID --Patternfile
com_group_for_gene = parent_parser.add_argument_group('Common options for gene information')
com_group_for_gene.add_argument("--repeatName", default=None, help="A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. 'all': all pre-defined genes will be analyzed;");
com_group_for_gene.add_argument("--UserDefinedUniqID", default=None, help="The name for storing results. Default: Pcr1. Must be different when simultaneously running RepeatHMM many times with the same setting for a same file");
com_group_for_gene.add_argument("--Patternfile", default=None, help="The file storing all predefined microsatellites (e.g., './reference_sts/hg38/hg38.trf.bed'). '.bed' for bed files and '.pa' for pa file. More than one files can be provided seperated by ';'. ");
com_group_for_gene.add_argument("--UserDefinedRepeat", default=UserDefinedRepeatdefault, help="The repeat information defined by users. If this option is given, the default repeat information will be revised. Default: ///////");

####################### --SplitAndReAlign --TRFOptions --minTailSize --minRepBWTSize --RepeatTime --BWAMEMOptions
com_group_for_splitalign = parent_parser.add_argument_group('Common options for re-alignment after splitting long reads using repeat regions')
com_group_for_splitalign.add_argument("--SplitAndReAlign", type=int, choices=[0, 1, 2], default=1, help="Split long reads using repeat region in long reads and re-align the non-repeat regions using BWA MEM. Default=0: not use SplitAndReAlign; 1: use SplitAndReAlign only; 2: combine 0 and 1")
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
com_group_for_hmm.add_argument("--emissionm", default=None, help="User-specified emission matrix for HMM. The number of rows must be the same as the 3*L+1 where L is the size of repeat unit and the number of columns must be 5. Probabilities in a row is separated by ',' and their sum must be 1. Probabilities of different rows are separated by ';'. Please pay more attention when providing this parameter. For CG repeat, the example of this matrix is '0.2,0.2,0.2,0.2,0.2;0.005,0.985,0.005,0.005,0;0.005,0.005,0.985,0.005,0;0.25,0.25,0.25,0.25,0;0.25,0.25,0.25,0.25,0;0.005,0.005,0.985,0.005,0;0.005,0.985,0.005,0.005,0'. Setting this option will override the setting of --hmm_insert_rate, --hmm_del_rate, --hmm_sub_rate and --SeqTech for emission matrix.");

#####################################
# --Onebamfile --SepbamfileTemp
parser_bam = subparsers.add_parser('BAMinput', parents=[parent_parser], help="Detect trinucleotide repeats from a BAM file", description="Detect trinucleotide repeats from a BAM file", epilog="For example, \n \
python %(prog)s --Onebamfile XXX.bam --repeatName HTT --hgfile XXX.fa \n \
python %(prog)s --Sepbamfile XXX%%s.bam --repeatName HTT --GapCorrection 1 FlankLength 30 --hgfile XXX.fa \n \
python %(prog)s --Onebamfile freeze4-all-merge.sort.bam --repeatName HTT --hgfile XXX.fa \n \
Final results was stored in logbam/RepBAM_*.log \n \
Note!!!!!! \n \
BAM file should be produced with the same version of 'hg' if '--hg' or '--hgfile' is specified. \n \
", formatter_class=RawTextHelpFormatter)

bamgroup = parser_bam.add_mutually_exclusive_group(); #required=True)
bamgroup.add_argument("--Onebamfile", default=None, help="A BAM file storing all alignments")
bamgroup.add_argument("--SepbamfileTemp", default=None, help="A separated BAM file template storing all alignments; separated by chromosome ids. For example, '--SepbamfileTemp'=mybam_Chr%%s_sorted.bam where '%%s' is chromosome id (1....22,x,y)")

parser_bam.set_defaults(func=BAMinput)

###################################### --fastq
parser_fasta = subparsers.add_parser('FASTQinput', parents=[parent_parser], help="Detect trinucleotide repeats from a FASTQ file", description="Detect trinucleotide repeats from a FASTQ file", epilog="For example, \n \
python %(prog)s --fastq XXX.fq --repeatName HTT --hgfile XXX.fa \n \
python %(prog)s --fastq XXX.fq --repeatName HTT --GapCorrection 1 --FlankLength 30 --hgfile XXX.fa \n \
python %(prog)s --fastq XXX.fq --repeatName HTT --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID PCR1 --hgfile XXX.fa \n \
python %(prog)s --repeatName atxn3 --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID sca3_pcr25_raw_test --fastq atxn3_data/rawdata/sam025.raw.fastq --hgfile XXX.fa \n \
Final results was stored in logfq/RepFQ_*.log \n \
", formatter_class=RawTextHelpFormatter)

parser_fasta.add_argument("--fastq", help="The file name for fasta sequences");
parser_fasta.set_defaults(func=FASTQinput)

############################################## --region 
parser_scan = subparsers.add_parser('Scan', parents=[parent_parser], help="For scanning whole data", description="Scan for a whole genome or a whole region given", epilog="For example, \n \
python %(prog)s --region chr4 --GapCorrection 1 --FlankLength 30 ; \n \
python %(prog)s --region chr4 --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID PCR1; \n \
python %(prog)s --region chr4 --GapCorrection 1 --FlankLength 30 \n \
python %(prog)s --region chr4 --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID pcr1 \n \
Final results was stored in logscan/scan_res/*.log \n \
", formatter_class=RawTextHelpFormatter)

parser_scan.add_argument("--region", default=None, help="The region where microsatellites are located. 'All' means all microsatellites (<10 nucleotides in repeat units), 'chrZ:X:Y' indicates a chromosome region where all microsatellites would be automatically searche. 'None': will only detect microsatellites predined in 'repeatName'");
parser_scan.add_argument("--thread", default=1, type=int, help="How many additional threads are used. Default: 1");
#parser_scan.add_argument("--outFolder", default='align/', help="How many additional threads are used. Default: 1");
parser_scan.add_argument("--conted", default=0, type=int, help="Whether continue the running last time. Default: 0");

bam2group = parser_scan.add_mutually_exclusive_group(); #required=True)
bam2group.add_argument("--Onebamfile", default=None, help="A BAM file storing all alignments")
bam2group.add_argument("--SepbamfileTemp", default=None, help="A separated BAM file template storing all alignments; separated by chromosome ids. For example, '--SepbamfileTemp'=mybam_Chr%%s_sorted.bam where '%%s' is chromosome id (1....22,x,y)")

parser_scan.set_defaults(func=scan)


if len(sys.argv)<2:
   parser.print_help();
else:
   args = parser.parse_args()
   args.func(args);





