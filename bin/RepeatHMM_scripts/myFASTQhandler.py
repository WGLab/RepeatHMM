
import os;
import sys;
import string;
import math;

import random;
import time
import resource
import numpy;

import logging

from . import myBAMhandler
from . import getAlignment
#from .myheader import *
from . import myheader
from . import myHMM
from . import printHMMmatrix
from . import myRepeatReAlignment

def get3part(mgloc, gene_start_end, repeat_start_end, repeatName, unique_file_id, analysis_file_id, hgfn, specifiedOptions): #
   logging.info('The region is %s from %d to %d' % (mgloc[0], gene_start_end[0], gene_start_end[1]))

   predata, mfasta, sufdata = myBAMhandler.getGene(repeatName, mgloc[0], gene_start_end, unique_file_id, analysis_file_id, hgfn, 10, specifiedOptions)

   upstreamstr = mfasta[:(repeat_start_end[0]-gene_start_end[0])]
   repregion   = mfasta[(repeat_start_end[0]-gene_start_end[0]):-(gene_start_end[1]-repeat_start_end[1])]
   downstreamstr = mfasta[-(gene_start_end[1]-repeat_start_end[1]):]

   return [upstreamstr, repregion, downstreamstr]


def getSCA3ForGivenGene(commonOptions, specifiedOptions, moreOptions):
	predres = []

	mgloc = moreOptions['mgloc']
	repeatName = moreOptions['repeatName']
	gene_start_end = moreOptions['gene_start_end']
	repeat_start_end = moreOptions['repeat_start_end']

	fastafile = specifiedOptions['fastafile']
 	unique_file_id = specifiedOptions['unique_file_id']
	analysis_file_id = specifiedOptions['analysis_file_id']

	hgfile = commonOptions['hgfile']
	MinSup = commonOptions['MinSup']

	repPat = moreOptions['repPat']

	myHMM.produce_for_repPat(commonOptions, moreOptions)
	len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)
	logging.info("len_repPat="+str(len_repPat))
	repPat = moreOptions['repPat']

	upstreamstr, repregion, downstreamstr = get3part(mgloc, gene_start_end, repeat_start_end, repeatName , unique_file_id, analysis_file_id, hgfile, specifiedOptions)

	if len(repregion)==0:
		logging.error("Not repeat region! please check!!"+repeatName+(' gene_location=[%d, %d], repeat_location=[%d, %d]' % (gene_start_end[0], gene_start_end[1], repeat_start_end[0], repeat_start_end[1])));
		sys.exit(1);

	logging.info("Test "+repeatName+(' gene_location=[%d, %d], repeat_location=[%d, %d]; upstreamsize=%d, downstreamsize=%d' % (gene_start_end[0], gene_start_end[1], repeat_start_end[0], repeat_start_end[1], repeat_start_end[0]-gene_start_end[0], gene_start_end[1]-repeat_start_end[1])))
	logging.info("Normal/Pathogenic repeats: %s" % mgloc[5])

	orirepeat = int(len(repregion)/float(len_repPat)) #3)

	logging.info("Orignal Test read="+'<<<'+repregion+'>>>'+(" #repeat=%d; #len=%d" % (orirepeat, len(repregion))))

	bwamem_w_option = 90*4;
	max_w_option, min_w_option = 500, 100;
	if bwamem_w_option<min_w_option: bwamem_w_option = min_w_option
	if bwamem_w_option>max_w_option: bwamem_w_option = max_w_option
	bwamem_w_option = bwamem_w_option + int(len(upstreamstr+repregion+downstreamstr)*0.4)
	if bwamem_w_option>max_w_option: bwamem_w_option = max_w_option

	start_time = time.time();

	bamfile = fastafile + '.bam'
	bamfile = fastafile + unique_file_id + '.bam'
	specifiedOptions['bamfile'] = bamfile

	myret = {}; myretdetail = {}

	#cmd = 'bwa mem -k17 -w'+str(bwamem_w_option)+' -W40 -r10 -A1 -B1 -O1 -E1 -L1 -t '+mthreads+' -v 2 '+hg_reference_and_index+'/'+hgfile+' '+ fastafile +' | samtools view -S -b | samtools sort > '+bamfile
	cmd = 'bwa mem -k17 -w'+str(bwamem_w_option)+' -W40 -r10 -A1 -B1 -O1 -E1 -L1 -t '+myheader.mthreads+' -v 2 '+hgfile+' '+ fastafile +' | samtools view -S -b | samtools sort > '+bamfile
	logging.info(cmd);
	os.system(cmd);
		
	cmd = 'samtools index '+bamfile
	logging.info(cmd);
	os.system(cmd)

	if (commonOptions['SplitAndReAlign'] in [0,2]) or testall:
		start_time = time.time();	
		if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2): print ('p2bamhmm start'); sys.stdout.flush()
		p2bamhmm = myBAMhandler.getRepeatForGivenGene(commonOptions, specifiedOptions, moreOptions)
		memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
		if p2bamhmm==None:
			if (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
				print ('ERROR None detection', moreOptions['repeatName'], mgloc)
				logging.error('ERROR None detection: ' + str( moreOptions['repeatName']) + ' ' + str(mgloc))
		
		myBAMhandler.addSumForAGene(p2bamhmm, myret, myretdetail, 'p2bamhmm', 2)
		end_time = time.time();
		if commonOptions['outlog'] <= myheader.M_WARNING and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2): print ('p2bamhmm end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()

	if ((commonOptions['SplitAndReAlign'] in [1,2]) or testall) and (commonOptions['SeqTech'] not in ["Illumina"]):
		start_time = time.time();
		if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2): print ('start p2sp'); sys.stdout.flush()

		#moreOptions['fafqfile'] = specifiedOptions['fastafile']
		#moreOptions['fafqtype'] = 'fq'
		moreOptions['fafqfile'] = bamfile
		moreOptions['fafqtype'] = 'bam'

		p2sp = myRepeatReAlignment.getRepeatCounts(commonOptions, specifiedOptions, moreOptions)
		memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
		if p2sp==None:
			if (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
				print ('ERROR None detection (sp)', moreOptions['repeatName'], mgloc)
				logging.error('ERROR None detection (sp): ' + str( moreOptions['repeatName']) + ' ' + str(mgloc))
		
		myBAMhandler.addSumForAGene(p2sp, myret, myretdetail, 'p2sp', 2)
		end_time = time.time();
		if commonOptions['outlog'] <= myheader.M_WARNING and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2): print ('p2sp end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()
	
	os.system('rm '+bamfile);
	os.system('rm '+bamfile+'.bai');

	return [myret, myretdetail];

def getSCA3forKnownGeneWithPartialRev(commonOptions, specifiedOptions):
	repeatName = commonOptions['repeatName'].lower()
	newinfo = commonOptions['specifiedRepeatInfo']

	moreOptions = myBAMhandler.get_gLoc(repeatName, commonOptions)
	moreOptions['repeatName'] = commonOptions['repeatName']

	res = getSCA3ForGivenGene(commonOptions, specifiedOptions, moreOptions);
	
	return res


if __name__=="__main__":
	
	LOG_FILENAME = "logsca3/mytest.txt"

	logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

	repPat = 'CTG';  forw_rerv = '+'
	repeatName = 'atxn3'; ra = [35, 37]  #  35, 57
	curreadfile = 'sim_data/atxn3_defpcr1_ins0.12_del0.02_sub0.02_cov100/1.fastq'
	gene_location = [92070888, 92072403]; repeat_location=[92071011, 92071052]
	unique_file_id, analysis_file_id = 'mytest_unique', 'mytest_analysis'
	
	curreadfile = 'sim_data/atxn3_ins0.12_del0.02_sub0.02_cov1900/1.fastq'; ra = [69, 70]
	curreadfile = 'sim_data/atxn3_ins0.12_del0.02_sub0.02_cov1600/55.fastq'; ra = [43, 81]



