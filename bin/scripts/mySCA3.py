
import os;
import sys;
import string;
import math;

import random;
import time
import resource
import numpy;

import logging

import trinucleotideRepeatRealSimulation
import findTrinucleotideRepeats
import getAlignment
from myheader import *
import myHMM
import useUnsymmetricalAlign
import printHMMmatrix
import myRepeatReAlignment


def getSCA3ForGivenGene(commonOptions, specifiedOptions, moreOptions):
	predres = []

	mgloc = moreOptions['mgloc']
	repeatgene = moreOptions['repeatgene']
	gene_start_end = moreOptions['gene_start_end']
	repeat_start_end = moreOptions['repeat_start_end']

	fastafile = specifiedOptions['fastafile']
 	unique_file_id = specifiedOptions['unique_file_id']
	analysis_file_id = specifiedOptions['analysis_file_id']
	isUnsymAlign = specifiedOptions['isUnsymAlign']

	isRemInDel = commonOptions['isRemInDel']
	isupdown = commonOptions['isupdown']
	isExtend = commonOptions['isExtend']
	hgfile = commonOptions['hgfile']
	MinSup = commonOptions['MinSup']

	repPat = moreOptions['repPat']

	myHMM.produce_for_repPat(commonOptions, moreOptions)
	len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)
	logging.info("len_repPat="+str(len_repPat))
	repPat = moreOptions['repPat']

	mgloc[1] = int(mgloc[1]);  mgloc[2] = int(mgloc[2]);
	mgloc[3] = int(mgloc[3]);  mgloc[4] = int(mgloc[4]);

	curgenestart =[mgloc[1], mgloc[2]]
	currepstart = [mgloc[3], mgloc[4]]

	upstreamstr, repregion, downstreamstr = trinucleotideRepeatRealSimulation.get3part(mgloc, gene_start_end, repeat_start_end, repeatgene , unique_file_id, analysis_file_id, hgfile)

	if len(repregion)==0:
		logging.error("Not repeat region! please check!!"+repeatgene+(' gene_location=[%d, %d], repeat_location=[%d, %d]' % (gene_start_end[0], gene_start_end[1], repeat_start_end[0], repeat_start_end[1])));
		sys.exit(1);

	logging.info("Test "+repeatgene+(' gene_location=[%d, %d], repeat_location=[%d, %d]; upstreamsize=%d, downstreamsize=%d' % (gene_start_end[0], gene_start_end[1], repeat_start_end[0], repeat_start_end[1], repeat_start_end[0]-gene_start_end[0], gene_start_end[1]-repeat_start_end[1])))
	logging.info("Normal/Pathogenic repeats: %s" % mgloc[7])

	orirepeat = int(len(repregion)/float(len_repPat)) #3)

	logging.info("Orignal Test read="+'<<<'+repregion+'>>>'+(" #repeat=%d; #len=%d" % (orirepeat, len(repregion))))

	bwamem_w_option = 90*4;
	max_w_option, min_w_option = 1000, 100;
	if bwamem_w_option<min_w_option: bwamem_w_option = min_w_option
	if bwamem_w_option>max_w_option: bwamem_w_option = max_w_option
	bwamem_w_option = bwamem_w_option + int(len(upstreamstr+repregion+downstreamstr)*0.4)
	if bwamem_w_option>max_w_option: bwamem_w_option = max_w_option

	start_time = time.time();

	moreOptions['gene_start_end'] = curgenestart
	moreOptions['repeat_orig_start_end'] = currepstart

	bamfile = fastafile + '.bam'
	bamfile = fastafile + unique_file_id + '.bam'
	specifiedOptions['bamfile'] = bamfile

	myret = {}; myretdetail = {}

	if testall or (not isUnsymAlign):
		start_time = time.time();
		print 'p2bamself start'; 
		cmd = 'bwa mem -k17 -w'+str(bwamem_w_option)+' -W40 -r10 -A1 -B1 -O1 -E1 -L1 -t 4 '+hg_reference_and_index+'/'+hgfile+' '+ fastafile +' | samtools view -S -b | samtools sort > '+bamfile
		logging.info(cmd);
		os.system(cmd);
		
		cmd = 'samtools index '+bamfile
		logging.info(cmd);
		os.system(cmd)
		
		p2bamself = findTrinucleotideRepeats.getRepeatForGivenGene2(commonOptions, specifiedOptions, moreOptions)
		memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
		findTrinucleotideRepeats.addSumForAGene(p2bamself, myret, myretdetail, 'p2bamself', 2)
		end_time = time.time();
		print ('p2bamself end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()

	if ((not isUnsymAlign) and (commonOptions['SplitAndReAlign'] in [0,2])) or testall:
		start_time = time.time();
		print 'p2bamhmm start'; sys.stdout.flush()
		p2bamhmm = findTrinucleotideRepeats.getRepeatForGivenGene(commonOptions, specifiedOptions, moreOptions)
		memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
		findTrinucleotideRepeats.addSumForAGene(p2bamhmm, myret, myretdetail, 'p2bamhmm', 2)
		end_time = time.time();
		print ('p2bamhmm end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()

	if testall or (not isUnsymAlign):
		os.system('rm '+bamfile);
		os.system('rm '+bamfile+'.bai');

	if ((not isUnsymAlign) and (commonOptions['SplitAndReAlign'] in [1,2])) or testall:
		start_time = time.time();
		print 'start p2sp'; sys.stdout.flush()
		moreOptions['fafqfile'] = specifiedOptions['fastafile']
		moreOptions['fafqtype'] = 'fq'
		p2sp = myRepeatReAlignment.getRepeatCounts(commonOptions, specifiedOptions, moreOptions)
		memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
		findTrinucleotideRepeats.addSumForAGene(p2sp, myret, myretdetail, 'p2sp', 2)
		end_time = time.time();
		print ('p2sp end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()

	if isUnsymAlign or testall:
		start_time = time.time();
		print 'start p2unsym'; sys.stdout.flush(); ###moreOptions['repPat']
		p2unsym = useUnsymmetricalAlign.useUnsymmetricalAlign(upstreamstr, repregion, downstreamstr, fastafile, repeatgene, moreOptions['repPat'], moreOptions['forw_rerv'], isRemInDel, isupdown, isExtend, bwamem_w_option, True, None, False, MinSup, commonOptions);
		memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
		findTrinucleotideRepeats.addSumForAGene(p2unsym, myret, myretdetail, 'p2unsym', 0)
		end_time = time.time();
		print ('p2unsym end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()

	return [myret, myretdetail];

def getSCA3forKnownGeneWithPartialRev(commonOptions, specifiedOptions):
	moreOptions = {}

	gLoc = commonOptions['gLoc']
	repeatgene = commonOptions['repeatgene'].lower()
	newinfo = commonOptions['specifiedGeneInfo']

	infospt = newinfo.split('/')
	repeatgene = repeatgene.lower()
	if gLoc.has_key(repeatgene):
		logging.info(repeatgene)
		mgloc = gLoc[repeatgene]
	else: mgloc = ['','','', '','', '','','']

	if len(infospt)<len(mgloc):
		logging.error("Error: wrong input for the gene of interes: %s whose length is %d less than the expected length %d" % (newinfo, len(infospt), len(mgloc)));
		sys.exit(101);

	for i in range(len(mgloc)):
		curnew = string.strip(infospt[i]);
		if not curnew=='': 
			mgloc[i] = curnew
	errorstr = ''
	for i in range(len(mgloc)):
		if string.strip(mgloc[i])=='' and i<len(mgloc)-1:
			errorstr += ('Error no information for %d\n' % i)
	if not errorstr=='':
		logging.error(errorstr);
		sys.exit(102);

	gene_start_end = [int(mgloc[1]), int(mgloc[2])]
	repeat_start_end = [int(mgloc[3]), int(mgloc[4])]

	moreOptions['gene_start_end'] = gene_start_end
	moreOptions['repeat_start_end'] = repeat_start_end
	moreOptions['mgloc'] = mgloc
	moreOptions['repeatgene'] = commonOptions['repeatgene']
	moreOptions['chr'] = mgloc[0]
	moreOptions['repPat'] = mgloc[5];
	moreOptions['forw_rerv']= mgloc[6]

	res = getSCA3ForGivenGene(commonOptions, specifiedOptions, moreOptions);
	
	return res


if __name__=="__main__":
	
	LOG_FILENAME = "logsca3/mytest_useUnsymmetricalAlign.txt"

	logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

	gLoc = findTrinucleotideRepeats.getDiseseGeneInRefGenomeLocation('hg38')
	
	repPat = 'CTG';  forw_rerv = '+'
	repeatgene = 'atxn3'; mgloc = gLoc[repeatgene]; ra = [35, 37]  #  35, 57
	curreadfile = 'sim_data/atxn3_defpcr1_ins0.12_del0.02_sub0.02_cov100/1.fastq'
	gene_location = [92070888, 92072403]; repeat_location=[92071011, 92071052]
	unique_file_id, analysis_file_id = 'mytest_useUnsymmetricalAlign_unique', 'mytest_useUnsymmetricalAlign_analysis'
	
	curreadfile = 'sim_data/atxn3_ins0.12_del0.02_sub0.02_cov1900/1.fastq'; ra = [69, 70]
	curreadfile = 'sim_data/atxn3_ins0.12_del0.02_sub0.02_cov1600/55.fastq'; ra = [43, 81]
	gene_location = [int(mgloc[1]), int(mgloc[2])]
	repeat_location = [int(mgloc[3]), int(mgloc[4])]
	
	gene_location = [int(mgloc[3])-2150, int(mgloc[4])+2150]
	
	upstreamstr, repregion, downstreamstr = get3part(mgloc, gene_location, repeat_location, repeatgene , unique_file_id, analysis_file_id)

	#print useUnsymmetricalAlign(upstreamstr, repregion, downstreamstr, curreadfile, repeatgene, repPat, forw_rerv, 0, 0, 3, ra, 400)
	#print 'dif'
	print useUnsymmetricalAlign(upstreamstr, repregion, downstreamstr, curreadfile, repeatgene, repPat, forw_rerv, 1, 18, 0, ra, 400)


