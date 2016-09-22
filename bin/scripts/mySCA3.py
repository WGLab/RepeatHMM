
import os;
import sys;
import string;
import math;

import random;
import time
import numpy;

import logging

import trinucleotideRepeatRealSimulation
import findTrinucleotideRepeats
import getAlignment
from UnsymmetricPairAlignment import UnsymmetricPairAlignment
from myheader import *
import myHMM


def useUnsymmetricalAlign(upstreamstr, repregion, downstreamstr, curreadfile, repeatgene, repPat, forw_rerv, isAlign, isupdown, isExtend, bandwo):
	pdebug = False; #pdebug=True; 
	mdebug = False; #mdebug=True;
	if mdebug:
		print 'upstreamstr', len(upstreamstr), upstreamstr
		print 'repregion', len(repregion), repregion
		print 'downstreamstr', len(downstreamstr), downstreamstr

	templatestr = upstreamstr+repregion+downstreamstr;

	bandw = bandwo; # + int(len(templatestr)*0.2)

	if mdebug or pdebug: print '\ntemp_all  :', len(templatestr[len(upstreamstr)-10:-len(downstreamstr)]), templatestr[len(upstreamstr)-10:-len(downstreamstr)], 'bandw=', bandw
	allsimreads = findTrinucleotideRepeats.myReadTxtFile(curreadfile)

	repregion_len_threhold = 3;
	repeatbeforeafter = isupdown
	toleratebeforeafter = 30+isupdown;
	if toleratebeforeafter>len(upstreamstr): toleratebeforeafter = len(upstreamstr)-1
	if toleratebeforeafter>len(downstreamstr): toleratebeforeafter = len(downstreamstr)-1

	logging.info("bwamem_w_option=%d, repeatbeforeafter=%d, toleratebeforeafter=%d\n" % (bandwo, repeatbeforeafter, toleratebeforeafter))

	repeatlength = []; replen = len(repregion)
	bp = getAlignment.getBasePair();
	
	isprint = 0; #1;	
	#for j in range(1, len(allsimreads), 4):
	for j in range(len(allsimreads)-3, 0, -4):
		if mdebug: print "cur line=", j 
		querystr = allsimreads[j];
		querystr = string.strip(querystr)

		#querystr = getAlignment.getComplementary3(bp, querystr)

		#alignres = UnsymmetricPairAlignment.unsymmetricPairWiseAlignment(templatestr[::-1], len(templatestr), querystr[::-1], len(querystr), match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandw, 0).split(';')
		alignres = UnsymmetricPairAlignment.unsymmetricPairWiseAlignment(templatestr, len(templatestr), querystr, len(querystr), match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandw, isprint).split(';')
		#alignres[0] = alignres[0][::-1]
		#alignres[1] = alignres[1][::-1]

		querystrR = getAlignment.getComplementary3(bp, querystr)
		alignresR = UnsymmetricPairAlignment.unsymmetricPairWiseAlignment(templatestr, len(templatestr), querystrR, len(querystrR), match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandw, isprint).split(';')

		#if len(alignresR[1])>len(alignres[1]):
		if int(alignresR[2])>int(alignres[2]):
			querystr = querystrR
			alignres = alignresR
		
		
		if mdebug: print 'done: '
		if mdebug or pdebug: print 'temp_first', len(templatestr), len(alignres[1]), len(alignres[1][len(upstreamstr)-10:-len(downstreamstr)]), alignres[1][len(upstreamstr)-10:-len(downstreamstr)]
		if mdebug or pdebug: print 'quey_first', len(querystr), len(alignres[0]), len(alignres[0][len(upstreamstr)-10:-len(downstreamstr)]), alignres[0][len(upstreamstr)-10:-len(downstreamstr)]
		res_temp = alignres[1].replace('-','');
		
		start_loc = templatestr.index(res_temp)
		if start_loc==-1:
			logging.error("Could not find the positio of aligned templates %s in %s" % res_temp, templatestr)
		
		if start_loc>=len(upstreamstr):
			logging.warning("Could not cover upstream (%d: %s) using aligned sequence %s for all %s" % (start_loc, upstreamstr, res_temp, templatestr))
		else:
			align_temp_i, query_i = 0, 0;
			start_repeat_loc = -1;
			end_repeat_loc = -1;
			temp_loc_start = [start_loc, 0, 0]; 
			consider3 = [upstreamstr, repregion, downstreamstr]
			beforenum = 0; afternum = 0;
			for i in range(len(temp_loc_start)):
				temp_i = temp_loc_start[i]
				curstreamstr = consider3[i]
				
				while temp_i<len(curstreamstr) and align_temp_i<len(alignres[1]):
					if alignres[1][align_temp_i]=='-':
						align_temp_i += 1;
						query_i += 1
					else:
						if not alignres[1][align_temp_i]==curstreamstr[temp_i]:
							logging.error("%d: template upstream (%s at %d) is different from aligned template (%s at %d)" % (i, curstreamstr[temp_i], temp_i, alignres[1][align_temp_i], align_temp_i))
						temp_i += 1
						align_temp_i += 1;
						query_i += 1
					if i==0 and temp_i<len(curstreamstr)-toleratebeforeafter:
						start_repeat_loc = query_i
					elif i==1 or (i==2 and temp_i<toleratebeforeafter):
						end_repeat_loc = query_i
					if align_temp_i<len(alignres[0]) and (not (alignres[0][align_temp_i]=='-') ):
						if i==0: beforenum += 1
						if i==2: afternum += 1
			if  beforenum<isupdown or afternum<isupdown:
                                logging.warning("Partial cover: orignal [start_repeat_loc, end_repeat_loc]=[%d, %d], [beforenum, afternum]=[%d, %d]< %d" % (start_repeat_loc, end_repeat_loc, beforenum, afternum, isupdown))
                                logging.warning("Partial cover: query="+alignres[0][start_repeat_loc:end_repeat_loc+1])
                                logging.warning("Partial cover: templ="+alignres[1][start_repeat_loc:end_repeat_loc+1])
			if beforenum<isupdown: 
				logging.warning("Partial cover: befor="+alignres[1][:start_repeat_loc])
				start_repeat_loc = -1
			if afternum<isupdown: 
				logging.warning("Partial cover: after="+alignres[1][end_repeat_loc+1:])
				end_repeat_loc = -1

			if end_repeat_loc<len(upstreamstr)+len(repregion)+1:
				logging.warning("Could not cover the whole repeat regions %d" % end_repeat_loc);
	
			if (start_repeat_loc==-1 or end_repeat_loc==-1) or end_repeat_loc-start_repeat_loc<repregion_len_threhold:
				logging.warning("Could not find the alignment %d result. The repeat region is %d, %d, " % (j, start_repeat_loc, end_repeat_loc));
			else:
				res_query = alignres[0].replace('-','');
				
				repeat_start_end = [start_repeat_loc, end_repeat_loc+1]
				repeat_start_end[0] -= isExtend; repeat_start_end[1] += isExtend;
				if isExtend>0 and repeat_start_end[0]<0: repeat_start_end[0]=0
				if isExtend>0 and repeat_start_end[1]>len(res_query)-1:
					repeat_start_end[1] = len(res_query)-1
			
	
				if mdebug or pdebug: print 'temp_rep', len(alignres[1][repeat_start_end[0]:repeat_start_end[1]]), alignres[1][repeat_start_end[0]:repeat_start_end[1]]
				if mdebug or pdebug: print 'quer_rep', len(alignres[0][repeat_start_end[0]:repeat_start_end[1]]), alignres[0][repeat_start_end[0]:repeat_start_end[1]]
	
				detectregion = alignres[0][repeat_start_end[0]:repeat_start_end[1]].replace('-','');
				if len(detectregion)<repregion_len_threhold:
					logging.warning("Not enough repeat region is detected")
					logging.warning('temp_rep: %d, %s' % (len(alignres[1][repeat_start_end[0]:repeat_start_end[1]]), alignres[1][repeat_start_end[0]:repeat_start_end[1]]))
					logging.warning('quer_rep: %d, %s' % (len(alignres[0][repeat_start_end[0]:repeat_start_end[1]]), alignres[0][repeat_start_end[0]:repeat_start_end[1]]))
					continue;
	
				if mdebug or pdebug: print len(detectregion), repeat_start_end, detectregion
				
				if isAlign:
					newstr, pre0, predstats = findTrinucleotideRepeats.getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, detectregion)

					#newstr, pre0, predstats = findTrinucleotideRepeats.getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, detectregion, match, mismatch, gap_in_perf, gap_in_read, gap_before_after)
				else:
					newstr, pre0, predstats = myHMM.hmmpred(detectregion, repPat, forw_rerv, repeatbeforeafter)
				repeatlength.append(len(newstr)/3)
				
				if repeatlength[-1]>70:
					logging.info("More repeat %3d: %s" % (repeatlength[-1], detectregion))
					logging.info("            %3d: %s" % (repeatlength[-1], predstats))
					logging.info("            %3d: %s" % (repeatlength[-1], newstr))
				
	p2, allocr = findTrinucleotideRepeats.get2Peaks(repeatlength)
	
	
	return [p2, allocr];

			
def getSCA3ForGivenGene(mgloc, fastafile, isUnsymAlign, unique_file_id, analysis_file_id, repeatgene, gene_start_end, repeat_start_end, isAlign=1, isupdown=90, isExtend=0):
	
	mgloc[1] = int(mgloc[1]);  mgloc[2] = int(mgloc[2]);
	mgloc[3] = int(mgloc[3]);  mgloc[4] = int(mgloc[4]);

	#curupdown = 120;   # curgenestart, currepstart
	#if curupdown<isupdown: curupdown=isupdown
	#curgenestart = [mgloc[3]-curupdown, mgloc[4]+curupdown]
	#if curgenestart[0]<1: curgenestart[0] = 1;
	#if curgenestart[0]<gene_start_end[0]: curgenestart[0]=gene_start_end[0]
	#if curgenestart[1]>gene_start_end[1]: curgenestart[1]=gene_start_end[1]
	curgenestart =[mgloc[1], mgloc[2]]
	
	currepstart = [mgloc[3], mgloc[4]]

	upstreamstr, repregion, downstreamstr = trinucleotideRepeatRealSimulation.get3part(mgloc, gene_start_end, repeat_start_end, repeatgene , unique_file_id, analysis_file_id)
	predres = []

	if len(repregion)==0:
		logging.error("Not repeat region! please check!!"+repeatgene+(' gene_location=[%d, %d], repeat_location=[%d, %d]' % (gene_start_end[0], gene_start_end[1], repeat_start_end[0], repeat_start_end[1])));
		sys.exit(1);

	logging.info("SCA3 "+repeatgene+(' gene_location=[%d, %d], repeat_location=[%d, %d]; upstreamsize=%d, downstreamsize=%d' % (gene_start_end[0], gene_start_end[1], repeat_start_end[0], repeat_start_end[1], repeat_start_end[0]-gene_start_end[0], gene_start_end[1]-repeat_start_end[1])))
	logging.info("Normal/Pathogenic repeats: %s" % mgloc[7])

	orirepeat = int(len(repregion)/3)

	logging.info("Orignal SCA3 read="+'<<<'+repregion+'>>>'+(" #repeat=%d; #len=%d" % (orirepeat, len(repregion))))

	bwamem_w_option = 90*4;
	max_w_option, min_w_option = 1000, 100;
	if bwamem_w_option<min_w_option: bwamem_w_option = min_w_option
	if bwamem_w_option>max_w_option: bwamem_w_option = max_w_option
	bwamem_w_option = bwamem_w_option + int(len(upstreamstr+repregion+downstreamstr)*0.4)

	start_time = time.time();

	bamfile = fastafile + '.bam'

	if True: #not isUnsymAlign:
		cmd = 'bwa mem -k17 -w'+str(bwamem_w_option)+' -W40 -r10 -A1 -B1 -O1 -E1 -L1 -t 4 '+hg38_reference_and_index+'/hg38.fa '+ fastafile +' | samtools view -S -b | samtools sort > '+bamfile
		logging.info(cmd);
		os.system(cmd);
			
		cmd = 'samtools index '+bamfile
		os.system(cmd)

	if not isUnsymAlign:
		p2all = findTrinucleotideRepeats.getRepeatForGivenGene(mgloc[0], repeatgene, curgenestart, currepstart, bamfile, mgloc[5], mgloc[6], isAlign, isupdown, isExtend, unique_file_id, analysis_file_id)
		p2 = p2all[2];
	else:
		p2all = useUnsymmetricalAlign(upstreamstr, repregion, downstreamstr, fastafile, repeatgene, mgloc[5], mgloc[6], isAlign, isupdown, isExtend, bwamem_w_option);
		p2 = p2all[0];
		pbwa = findTrinucleotideRepeats.getRepeatForGivenGene2(mgloc[0], repeatgene, curgenestart, currepstart, bamfile, mgloc[5], mgloc[6], isAlign, isupdown, isExtend, unique_file_id, analysis_file_id)
		pbwa = pbwa[2]; pbwa.sort();
		if len(pbwa)==1: pbwa = [pbwa[0], pbwa[0]]
		if len(pbwa)==0: pbwa = [0, 0]
	
	p2.sort();
	if len(p2)==1: p2 = [p2[0], p2[0]];
	if len(p2)==0: p2 = [0, 0];

	if isUnsymAlign:
		logging.info("Result : %s, %s; Elapsed time=%.2f\n" % (str(p2[0]), str(p2[1]), (time.time() - start_time)))
	else:
		logging.info("Result : %s, %s; wrong alignment %d for %s. Elapsed time=%.2f\n" % (str(p2[0]), str(p2[1]), p2all[5], bamfile, (time.time() - start_time)))
	
	if not isUnsymAlign:
		predres.append([p2, p2all[5]])
	else: predres.append([p2, pbwa])

	return predres

def getSCA3forKnownGeneWithPartialRev(gLoc, fastafile, isUnsymAlign, unique_file_id, analysis_file_id, repeatgene, newinfo, isAlign=1, isupdown=90, isExtend=0):
        infospt = newinfo.split('/')
	repeatgene = repeatgene.lower()
	if gLoc.has_key(repeatgene):
		print gLoc.keys()
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
		if string.strip(mgloc[i])=='':
			errorstr += ('Error no information for %d\n' % i)
	if not errorstr=='':
		logging.error(errorstr);
		sys.exit(102);

	gene_start_end = [int(mgloc[1]), int(mgloc[2])]
	repeat_start_end = [int(mgloc[3]), int(mgloc[4])]

	#print mgloc; sys.exit(103);

        res = getSCA3ForGivenGene(mgloc, fastafile, isUnsymAlign, unique_file_id, analysis_file_id, repeatgene, gene_start_end, repeat_start_end, isAlign, isupdown, isExtend)
	
	return res


#                  0         1             2            3            4         5     6     7
#gloc['fmr1'] = ['chrX', '147911951', '147951127', '14799051', '147912110', 'CGG', '+', '6-53:230+/55-200']

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


