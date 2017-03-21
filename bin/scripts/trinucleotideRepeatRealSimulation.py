
import os;
import sys;
import string;
import math;

import random;
import time
import resource
import numpy;

import logging


import findTrinucleotideRepeats
import getAlignment
from myheader import *
import myHMM
import useUnsymmetricalAlign
import printHMMmatrix
import myRepeatReAlignment

def readList(fname, mtype=int):
	f = open(fname,'r')
	mlines = f.readlines();
	data = [];
	for ml in mlines:
		if ml[0]=='#': continue;
		mlsp = ml.split();
		curd = []
		for curi in mlsp:
			curi = string.strip(curi);
			if not curi=='': curd.append(mtype(curi));
		data.append(curd);
	f.close();
	return data;

def writeList(fname, mlist):
	f = open(fname, 'w')
	for ml in mlist:
		for ml1_ind in range(len(ml)):
			if ml1_ind==len(ml)-1: f.write(str(ml[ml1_ind])+'\n')
			else: f.write(str(ml[ml1_ind])+'\t')
	f.close();


def create_quality_for_pacbio_read(length, mean, sd):
    qualities = [fastq_Sanger_quality[min(max(0, int(random.gauss(mean, sd))), 93)] for i in range(length)]
    return ''.join(qualities)

def get3part(mgloc, gene_start_end, repeat_start_end, repeatgene, unique_file_id, analysis_file_id, hgfn): #
	logging.info('The region is %s from %d to %d' % (mgloc[0], gene_start_end[0], gene_start_end[1]))
	
	mfasta = findTrinucleotideRepeats.getGene(repeatgene, mgloc[0], gene_start_end, unique_file_id, analysis_file_id, hgfn)
	
	upstreamstr = mfasta[:(repeat_start_end[0]-gene_start_end[0])]
	repregion   = mfasta[(repeat_start_end[0]-gene_start_end[0]):-(gene_start_end[1]-repeat_start_end[1])]
	downstreamstr = mfasta[-(gene_start_end[1]-repeat_start_end[1]):]

	return [upstreamstr, repregion, downstreamstr]

def getMidPos(seq, na3):
	findpos = False;

	for pos in range(int(len(seq)*0.45), len(seq)-2):
		if seq[pos:(pos+len(na3))]==na3: #
			findpos = True;
			break;
	if findpos: return pos;
	else: return 0;

def getNA3fromRegion(seq, na3):
	patlen = len(na3);
	
	poslist = []; lenlist = []
	loci = 0; preloci = 0;
	while loci<len(seq)-patlen: #3:
		if seq[loci:loci+patlen] == na3:
			if preloci<loci:
				poslist.append(preloci); lenlist.append(loci-preloci);
			poslist.append(loci); lenlist.append(patlen);
			loci = loci + patlen;
			preloci = loci;
		else: loci += 1;
	if preloci<loci:
		poslist.append(preloci); lenlist.append(len(seq)-preloci);
	elif preloci<len(seq) and loci<len(seq):
		poslist.append(preloci); lenlist.append(len(seq)-preloci);

	plen = 0; #
	while True:
		while lenlist[plen]>patlen:
			poslist.insert(plen+1, poslist[plen]+patlen);
			lenlist.insert(plen+1, lenlist[plen]-patlen);
			lenlist[plen] = patlen;
			plen += 1;
		
		plen += 1;
		if plen >= len(poslist): break;
	
	return [poslist, lenlist];

def getRandomNA(bp):
	#print bp;
	prob = int(random.uniform(1,len(bp)))-1
	return bp[prob];

def mutateStr(seq, insert_rate, del_rate, sub_rate, bp5):
	seqlen = len(seq);
	insert_num = int(seqlen*insert_rate); #if insert_num<1: insert_num = 1;
	del_num = int(seqlen*del_rate); #if del_num<1: del_num = 1;
	sub_num = int(seqlen*sub_rate); #if sub_num<1: sub_num = 1;

	mmin = 0;
	insert_num = int(random.gauss(insert_num, 2)+0.5); 
	if insert_num<mmin: insert_num = mmin;
	del_num = int(random.gauss(del_num, 1)+0.5); 
	if del_num<mmin: del_num = mmin;
	sub_num = int(random.gauss(sub_num, 1)+0.5); 
	if sub_num<mmin: sub_num = mmin;

	#insertion could occur in the same position
	insert_pos = [random.randint(0,seqlen-1) for _ in xrange(insert_num)]; insert_pos.sort();
	allpos = range(seqlen);
	if del_num==0: del_pos=[] # deletion could not occur in the same position
	else: del_pos = random.sample(allpos, del_num); del_pos.sort();
	if sub_num==0: sub_pos = [] # substitution could occur in the smae position
	else: sub_pos = random.sample(allpos, sub_num); sub_pos.sort();

	newstr = ''
	curp = 0; ij, dj, sj = 0, 0, 0;
	while curp < seqlen:
		while ij<insert_num and insert_pos[ij]==curp:
			newstr += getRandomNA(bp5[''])
			ij += 1;
		
		if dj<del_num and del_pos[dj]==curp:
			dj += 1;
		else:
			while sj<sub_num and sub_pos[sj]<curp: sj += 1;

			if sj<sub_num and sub_pos[sj]==curp:
				newstr += getRandomNA(bp5[seq[curp]]);
				sj += 1;
			else: newstr += seq[curp]

		curp += 1;

	if ij<insert_num-1: logging.warning("Not all insert positions are used: used=%d, total=%d with largest position=%d and current position=%d for seqlen=%d" % (ij, insert_num, insert_pos[-1], insert_pos[ij], seqlen));
	if dj<del_num-1: logging.warning("Not all deletion positions are used: used=%d, total=%d with largest position=%d and current position=%d for seqlen=%d" % (dj, del_num, del_pos[-1], del_pos[dj], seqlen))
	if sj<sub_num-1: logging.warning("Not all substitution positions are used: used=%d, total=%d with largest position=%d and current position=%d for seqlen=%d" % (sj, sub_num, sub_pos[-1], sub_pos[sj], seqlen))

	return newstr;


def getMynormNum(avg, std,  mmin, mmax, mstr):
	currn = int(numpy.random.normal(avg, std, 1)[0])
	trytimes = 0;
	while (currn<mmin or currn>mmax):
		currn = int(numpy.random.normal(avg, std, 1)[0])
		trytimes = trytimes + 1
		if trytimes>1000:
			logging.error("Cannot find a number for "+mstr+(": trytimes=%d" % (trytimes)))
			logging.error(asd)
	return currn

def getrandomstr(newup, currepeat, newdown, shortennum=1):
	newrandstr = newup + currepeat + newdown
	lenall = len(newrandstr);
	lenfirst = len(newup)-shortennum; avg1 = lenfirst/2
	lensecond = len(newdown)-shortennum; avg2 = lensecond/2

	pos1 = getMynormNum(avg1, 10, 1, lenfirst-1, 'First random position')
	pos2 = getMynormNum(avg2, 10, 1, lensecond-1, 'Second random position')
	
	return newrandstr[pos1:(-pos2)]
def getNormalDist(avgnormrep, mmin, mmax):
	currn = getMynormNum(avgnormrep, 5, mmin, mmax, 'repeats')
	return currn

def getDefinedSizeDif(mmin, mmax, m_repeatSizeDif):
        if len(m_repeatSizeDif)==1:
                repeatSizeDif = [m_repeatSizeDif[0], m_repeatSizeDif[0]]
        else:
                repeatSizeDif = [m_repeatSizeDif[0], m_repeatSizeDif[1]]

        repeat1 = int(random.uniform(mmin+repeatSizeDif[1], mmax-repeatSizeDif[1]))
        repeat2array = []
        for curr2 in range(repeat1-repeatSizeDif[1], repeat1+repeatSizeDif[1]+1):
                if not (curr2>repeat1-repeatSizeDif[0] and curr2<repeat1+repeatSizeDif[0]):
                        repeat2array.append(curr2);
        rep2ind = random.uniform(0,len(repeat2array)-1)
        repeat2 = repeat2array[int(rep2ind+0.5)]

        retrep = [repeat1, repeat2]; retrep.sort();

        return retrep


def getSimForGivenGene(commonOptions, specifiedOptions, moreOptions):
	mgloc = moreOptions['mgloc']
	repeatgene = moreOptions['repeatgene']
	gene_start_end = moreOptions['gene_start_end']
	repeat_start_end = moreOptions['repeat_start_end']

	repeatrange = moreOptions['repeatrange']
	isPartial = moreOptions['isPartial']

	isRemInDel = commonOptions['isRemInDel']
	isupdown = commonOptions['isupdown']
	isExtend = commonOptions['isExtend']
	hgfile = commonOptions['hgfile']
	MinSup = commonOptions['MinSup']

	unique_file_id = specifiedOptions['unique_file_id']
	analysis_file_id = specifiedOptions['analysis_file_id']
	simulation_file_id = specifiedOptions['simulation_file_id']
	isUnsymAlign = specifiedOptions['isUnsymAlign']
	
	insert_rate = specifiedOptions['insert_rate']
	del_rate = specifiedOptions['del_rate']
	sub_rate = specifiedOptions['sub_rate']
	coverage = specifiedOptions['coverage']
	randTimes = specifiedOptions['randTimes']
	repeatSizeDif = specifiedOptions['repeatSizeDif']

	simfolder = 'sim_data/'
	if not os.path.isdir(simfolder): 
		os.system("mkdir "+simfolder);

	#simfile = simfolder + repeatgene + simulation_file_id
	simfile = simfolder + simulation_file_id
	if not os.path.isdir(simfile):
		os.system("mkdir "+simfile)

	avgnormrep = int(mgloc[6][1:])
	
	mgloc[1] = int(mgloc[1]);  mgloc[2] = int(mgloc[2]);
	mgloc[3] = int(mgloc[3]);  mgloc[4] = int(mgloc[4]);
	#curupdown = 120 #
	#if curupdown<isupdown: curupdown=isupdown
	#curgenestart = [mgloc[3]-curupdown, mgloc[4]+curupdown]
	#if curgenestart[0]<1: curgenestart[0] = 1;
	#if curgenestart[0]<gene_start_end[0]: curgenestart[0]=gene_start_end[0]
	#if curgenestart[1]>gene_start_end[1]: curgenestart[1]=gene_start_end[1]
	#currepstart = [mgloc[3], mgloc[4]]

	upstreamstr, repregion, downstreamstr = get3part(mgloc, gene_start_end, repeat_start_end, repeatgene , unique_file_id, analysis_file_id, hgfile)
	predres = []

	if len(repregion)==0:
		logging.error("Not repeat region! please check!!"+repeatgene+(' gene_location=[%d, %d], repeat_location=[%d, %d]' % (gene_start_end[0], gene_start_end[1], repeat_start_end[0], repeat_start_end[1])));
		sys.exit(1);

	na3 = getAlignment.getPattern(mgloc[5], mgloc[6])
	
	pos = getMidPos(repregion, na3);

	logging.info("Simulation for "+repeatgene+(' gene_location=[%d, %d], repeat_location=[%d, %d]; upstreamsize=%d, downstreamsize=%d' % (gene_start_end[0], gene_start_end[1], repeat_start_end[0], repeat_start_end[1], repeat_start_end[0]-gene_start_end[0], gene_start_end[1]-repeat_start_end[1])))
	logging.info("Normal/Pathogenic repeats: %s; avgnormrep=%d" % (mgloc[7], avgnormrep))
	logging.info("Detect a point in the middle for more repeat insertion");	
	logging.info(repregion+' '+repregion[pos:(pos+3)]+' '+str(pos)+' '+na3)
	logging.info(("Mutation info: insertion rate=%.2f, deletion rate=%.2f, substitution rate=%.2f, normal_range=[%d, %d]/path_range=[%d, %d]" % (insert_rate, del_rate, sub_rate, repeatrange[0][0], repeatrange[0][1], repeatrange[1][0], repeatrange[1][1])))

	bp = getAlignment.getBasePair(); bpkeys = bp.keys(); bpkeys.sort();
	bp5 = {}
	bp5[''] = bpkeys
	bp5['A'] = ['C', 'G', 'T']
	bp5['C'] = ['A', 'G', 'T']
	bp5['G'] = ['A', 'C', 'T']
	bp5['T'] = ['A', 'C', 'G']

	poslist, lenlist = getNA3fromRegion(repregion, na3);
	if True:
		posstr = ''
		for posind in range(len(poslist)):
			posstr += '|'
			for j in range(lenlist[posind]-1):
				posstr += ' '
	posind = range(1, len(poslist))
	
	orirepeat = int(len(repregion)/float(len(na3))) #3)

	logging.info("Orignal simulation read="+'<<<'+repregion+'>>>'+(" #repeat=%d; #len=%d" % (orirepeat, len(repregion))))

	bwamem_w_option = repeatrange[1][1]*4;
	max_w_option, min_w_option = 1000, 100;
	if bwamem_w_option<min_w_option: bwamem_w_option = min_w_option
	if bwamem_w_option>max_w_option: bwamem_w_option = max_w_option
	bwamem_w_option = bwamem_w_option + int(len(upstreamstr+repregion+downstreamstr)*0.5)
	if bwamem_w_option>max_w_option: bwamem_w_option = max_w_option

	logging.info("bwamem_w_option="+str(bwamem_w_option))

	previous_sim_repeats = []
	produced_repeat_file = simfile+'_produced_repeats.txt'
	if os.path.isfile(produced_repeat_file):
		previous_sim_repeats = readList(produced_repeat_file)
	cur_sim_repeats = []
	
	alwaysproduced = False; 
	alwaysproduced = True;
	randomproduced = False; #True; #False; #True;
	refcov = 10000; #9400
	refcov = 3000; 
	curcovstr = ('cov%d' % coverage)
	curcovind = simfile.index(curcovstr)
	refcovfolder = simfile[:curcovind]+('cov%d' % refcov)+simfile[(curcovind+len(curcovstr)):]
	if os.path.isfile(refcovfolder+'_produced_repeats.txt'):
		refcovrepeats = readList(refcovfolder+'_produced_repeats.txt')
	else: randomproduced = True;

	alwaysproduced = True;
	randomproduced = True;

	moreOptions['repeat_orig_start_end'] = [moreOptions['repeat_start_end'][0], moreOptions['repeat_start_end'][1]]
	moreOptions['chr'] = mgloc[0]
	moreOptions['repPat'] = mgloc[5];
	moreOptions['forw_rerv']= mgloc[6]

	myHMM.produce_for_repPat(commonOptions, moreOptions)
	len_repPat = printHMMmatrix.get_len_repPat(moreOptions['repPat'], commonOptions)
	logging.info("len_repPat="+str(len_repPat))

	mysimsum = {}
	
	for rt in range(randTimes):
		myret = {}; myretdetail = {}

		start_time = time.time();

		curreadfile = simfile + '/' + str(rt) +'.fastq'
		bamfile = curreadfile + '.bam'
		bamfile = curreadfile + unique_file_id + '.bam'
		
		specifiedOptions['bamfile'] = bamfile
		specifiedOptions['fastafile'] = curreadfile

		if isUnsymAlign and os.path.isfile(curreadfile) and rt<len(previous_sim_repeats) and not (randomproduced or alwaysproduced):
			currep2 = previous_sim_repeats[rt]
			cur_sim_repeats.append(currep2)
			
			logging.info("Simutlate repeats are %d and %d at %d round" % (currep2[0], currep2[1], rt))
			logging.info("File exist "+curreadfile)
		elif (not isUnsymAlign) and os.path.isfile(bamfile) and os.path.isfile(bamfile+'.bai') and rt<len(previous_sim_repeats) and rt<len(previous_sim_repeats):
			currep2 = previous_sim_repeats[rt]
			
			cur_sim_repeats.append(currep2)
			
			logging.info("Simutlate repeats are %d and %d at %d round" % (currep2[0], currep2[1], rt))
			logging.info("File exist "+bamfile)
		else: 
			if (randomproduced or alwaysproduced): #:
				if len(repeatSizeDif)==0:
					repeat1 = getNormalDist(avgnormrep, repeatrange[0][0], repeatrange[0][1])
					repeat1 = int(random.uniform(repeatrange[0][0], repeatrange[0][1]))
					repeat2 = int(random.uniform(repeatrange[1][0], repeatrange[1][1]))
				else:
					repeat1, repeat2 = getDefinedSizeDif(repeatrange[0][0], repeatrange[1][1], repeatSizeDif)
					
				currep2 = [repeat1, repeat2]; currep2.sort();
				cur_sim_repeats.append(currep2)
				
				logging.info("Simutlate repeats are %d and %d at %d round" % (currep2[0], currep2[1], rt))
			else: 
				currep2 = refcovrepeats[rt]
				cur_sim_repeats.append(currep2)
				
				refsq = findTrinucleotideRepeats.myReadTxtFile(refcovfolder+ '/' + str(rt) +'.fastq')
				
			
			allsimreads = [];
			currep2.sort(); num_long_reads = int(coverage/(1.5**(currep2[1]/float(currep2[0])-1))+0.5)
			logging.info("num_long_reads="+str(num_long_reads))

			shortennum = [1, (currep2[1]-currep2[0])*len(na3)/2]
			coverage2 = [coverage*currep2[1]/(currep2[1]+currep2[0]), coverage*currep2[0]/(currep2[1]+currep2[0])]
			logging.info("coverage2=[%d, %d] (total=%d) at %d round, shortennum=[%d, %d] for %d/%d" % (coverage2[0], coverage2[1], coverage, rt, shortennum[0], shortennum[1], currep2[0], currep2[1]))
			print ("coverage2=[%d, %d] (total=%d) at %d round, shortennum=[%d, %d] for %d/%d" % (coverage2[0], coverage2[1], coverage, rt, shortennum[0], shortennum[1], currep2[0], currep2[1]))
			for i_ind in range(len(currep2)):
				i = currep2[i_ind]
				for j in range(coverage2[i_ind]):

					if isPartial and i==currep2[1] and j>=num_long_reads: continue
					if not (randomproduced or alwaysproduced):
						if coverage*5<len(refsq):
							allsimreads.append(refsq[coverage*4]);
							allsimreads.append(refsq[coverage*4+1]);
							allsimreads.append(refsq[coverage*4+2]);
							allsimreads.append(refsq[coverage*4+3]);
							continue;

					if i < orirepeat:
						neworder = random.sample(posind, i);
						neworder.sort();
						newrepeat = ''
						for noi in neworder:
							newrepeat += repregion[poslist[noi]:(poslist[noi]+lenlist[noi])]
					else:
						curind = int(random.uniform(1, len(poslist)));
						if curind==0:
							newrepeat = (na3*(i-orirepeat)) + repregion
						elif curind==len(poslist)-1:
							newrepeat = repregion + (na3*(i-orirepeat))
						else:
							curpos = poslist[curind]
							newrepeat = repregion[:curpos] + (na3*(i-orirepeat)) + repregion[curpos:]
					
					newup = mutateStr(upstreamstr, insert_rate, del_rate, sub_rate, bp5)
					currepeat = mutateStr(newrepeat, insert_rate, del_rate, sub_rate, bp5)
					newdown =  mutateStr(downstreamstr, insert_rate, del_rate, sub_rate, bp5)
					
					if j<3:
						logging.debug('Perf: '+('%.0f' % (len(newrepeat)/3))+'>>>'+newrepeat)
						logging.debug('Mutat:'+('%.0f' % (len(currepeat)/3))+'>>>'+currepeat)
					
					allsimreads.append('@'+repeatgene+':'+str(i)+':'+str(j))
					if isPartial: #
						allsimreads.append(newup+currepeat+newdown)
					else:
						newrandstr = getrandomstr(newup, currepeat, newdown, shortennum[i_ind])
						allsimreads.append(newrandstr)
					allsimreads.append('+')
					allsimreads.append(create_quality_for_pacbio_read(len(allsimreads[-2]), 30, 10));
		
			findTrinucleotideRepeats.myWriteTxtFile(allsimreads, curreadfile)

		if testall or (not isUnsymAlign):
			start_time = time.time();
			print 'start p2bamself'; sys.stdout.flush()
			cmd = 'bwa mem -k17 -w'+str(bwamem_w_option)+' -W40 -r10 -A1 -B1 -O1 -E1 -L1 -t 4 '+hg_reference_and_index+'/'+hgfile+' '+ curreadfile +' | samtools view -S -b | samtools sort > '+bamfile
			logging.info(cmd); os.system(cmd);
			cmd = 'samtools index '+bamfile; os.system(cmd)

			p2bamself = findTrinucleotideRepeats.getRepeatForGivenGene2(commonOptions, specifiedOptions, moreOptions)
			memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
			findTrinucleotideRepeats.addSumForAGene(p2bamself, myret, myretdetail, 'p2bamself', 2)
			end_time = time.time();
			print ('p2bamself end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()

		if ((not isUnsymAlign) and (commonOptions['SplitAndReAlign'] in [0,2])) or testall:
			start_time = time.time();
			print 'start p2bamhmm'; sys.stdout.flush()
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
			print 'start p2unsym'; sys.stdout.flush()  ###
			p2unsym = useUnsymmetricalAlign.useUnsymmetricalAlign(upstreamstr, repregion, downstreamstr, curreadfile, repeatgene, moreOptions['repPat'], moreOptions['forw_rerv'], isRemInDel, isupdown, isExtend, bwamem_w_option, True, None, False, MinSup, commonOptions);
			memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
			findTrinucleotideRepeats.addSumForAGene(p2unsym, myret, myretdetail, 'p2unsym', 0)
			end_time = time.time();
			print ('p2unsym end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()

		mysimsum[rt] = [myret, myretdetail]

		if isUnsymAlign:
			logging.info("Result for simutlation repeats are %d and %d at %d round: %s, %s; Elapsed time=%.2f\n" % (currep2[0], currep2[1], rt, str(p2unsym[0]), str(p2unsym[1]), (time.time() - start_time)))

	if len(previous_sim_repeats)<=len(cur_sim_repeats):
		writeList(produced_repeat_file, cur_sim_repeats)

	return mysimsum

def getRange1(str1):
	strsp = str1.split('/');
	normrange = [100000,0];
	for sp1 in strsp:
		curr = sp1.split('-')
		if len(curr)==1:
			curr = curr[0]
			if curr[-1]=='+': curr = curr[:-1];
			curmin = int(curr);
			if curmin<normrange[0]: normrange[0] = curmin;
			if normrange[1]<curmin: normrange[1] = int(curmin * 1.25)
		elif len(curr)==2:
			curmin = int(curr[0]); curmax = int(curr[1]);
			if curmin<normrange[0]: normrange[0] = curmin;
			if curmax>normrange[1]: normrange[1] = curmax
	return normrange

def getRepeatRange2(mstr):
	strsp = mstr.split(':')
	normstr = strsp[0];     pathstr = strsp[1];
	normrange = getRange1(normstr);
	pathrange = getRange1(pathstr);
	return [normrange, pathrange]

def getRepeatRange(mstr):
	strsp = mstr.split(':')
	normstr = strsp[0];	pathstr = strsp[1];
	
	normrange = getRange1(normstr);
	pathrange = getRange1(pathstr);

	allrange = [normrange[0], pathrange[1]];
	if allrange[0]>pathrange[0]:
		logging.error("Too less repeat for Pathogenic:"+('Pathogenic=[%d, %d], Normal=[%d, %d]' % (pathrange[0], pathrange[1], normrange[0], normrange[1])))
	if pathrange[1]<normrange[1]: 
		logging.error("Too more repeat for Normal:"+('Pathogenic=[%d, %d], Normal=[%d, %d]' % (pathrange[0], pathrange[1], normrange[0], normrange[1])))

	if allrange[0]>=allrange[1]:
		logging.error("Wrong repeat range="+('[%d, %d]' % (allrange[0], allrange[1])))

	return allrange;

def getSimforKnownGeneWithPartialRev(commonOptions, specifiedOptions):
	moreOptions = {}

	gLoc = commonOptions['gLoc']
	repeatgene = commonOptions['repeatgene']
	newinfo = commonOptions['specifiedGeneInfo']

	infospt = newinfo.split('/')
	repeatgene = repeatgene.lower()
	if gLoc.has_key(repeatgene):
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
		if string.strip(mgloc[i])=='':
			errorstr += ('Error no information for %d\n' % i)
	if not errorstr=='':
		logging.error(errorstr);
		sys.exit(102);

	gene_start_end = [int(mgloc[1]), int(mgloc[2])]
	repeat_start_end = [int(mgloc[3]), int(mgloc[4])]
	repeatrange = getRepeatRange(mgloc[7])
	repeatrange = getRepeatRange2(mgloc[7])

	moreOptions['gene_start_end'] = gene_start_end
	moreOptions['repeat_start_end'] = repeat_start_end
	moreOptions['mgloc'] = mgloc
	moreOptions['repeatgene'] = commonOptions['repeatgene']
	moreOptions['repeatrange'] = repeatrange
	moreOptions['isPartial'] = True;

	res = getSimForGivenGene(commonOptions, specifiedOptions, moreOptions)
	
	return res

def getSimforKnownGene(commonOptions, specifiedOptions):
#                  0         1             2            3            4         5     6     7
#gloc['fmr1'] = ['chrX', '147911951', '147951127', '14799051', '147912110', 'CGG', '+', '6-53:230+/55-200']
	moreOptions = {}

	gLoc = commonOptions['gLoc']
	repeatgene = commonOptions['repeatgene']
	
	repeatgene = repeatgene.lower()
	mgloc = findTrinucleotideRepeats.get_gLoc(repeatgene, gLoc);

	repeatrange = getRepeatRange(mgloc[7])
	defaultupdownstreamsize = repeatrange[1]*25;
	maxupdown = 6000; minupdown = 4000;
	maxupdown = 1500; minupdown = 1500;
	if defaultupdownstreamsize>maxupdown: defaultupdownstreamsize = maxupdown
	elif defaultupdownstreamsize<minupdown: defaultupdownstreamsize = minupdown;

	gene_start_end = [int(mgloc[1]), int(mgloc[2])]
	repeat_start_end = [int(mgloc[3]), int(mgloc[4])]

	upstreampos = repeat_start_end[0] - defaultupdownstreamsize
	if upstreampos<1: upstreampos = 1
	downstreampos = repeat_start_end[1] + defaultupdownstreamsize

	gene_start_end = [upstreampos, downstreampos]
	repeatrange = getRepeatRange2(mgloc[7])

	moreOptions['gene_start_end'] = gene_start_end
	moreOptions['repeat_start_end'] = repeat_start_end
	moreOptions['mgloc'] = mgloc
	moreOptions['repeatgene'] = commonOptions['repeatgene']
	moreOptions['repeatrange'] = repeatrange
	moreOptions['isPartial'] = False;

	res = getSimForGivenGene(commonOptions, specifiedOptions, moreOptions)
	return res;


def getSim(commonOptions, specifiedOptions):
	res = []

	glkeys = gLoc.keys(); glkeys.sort()
	for glk in glkeys:
		commonOptions['repeatgene'] = glk
		res.extend(getSimforKnownGene(commonOptions, specifiedOptions))

	return res;

if __name__=="__main__":
	
	LOG_FILENAME = "log/mytest_useUnsymmetricalAlign.txt"

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

	print useUnsymmetricalAlign(upstreamstr, repregion, downstreamstr, curreadfile, repeatgene, repPat, forw_rerv, 1, 18, 0, ra, 400)


