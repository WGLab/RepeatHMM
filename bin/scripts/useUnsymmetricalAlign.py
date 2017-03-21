
import os;
import sys;
import string
import logging

import findTrinucleotideRepeats
import getAlignment
from UnsymmetricPairAlignment import UnsymmetricPairAlignment
import myHMM
from myheader import *
import printHMMmatrix

def useUnsymmetricalAlign(upstreamstr, repregion, downstreamstr, curreadfile, repeatgene, repPat, forw_rerv, isRemInDel, isupdown, isExtend, bandwo, seqbp, currep2, printTwoTail, MinSup, commonOptions):
	mdebug = False; #mdebug=True;
	if mdebug:
		print 'upstreamstr', len(upstreamstr), upstreamstr
		print 'repregion', len(repregion), repregion
		print 'downstreamstr', len(downstreamstr), downstreamstr

	templatestr = upstreamstr+repregion+downstreamstr;

	len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)
	logging.info("len_repPat="+str(len_repPat))

	print 'rep patterns', repPat, commonOptions['CompRep']

	bandw = bandwo; 

	if mdebug: print '\ntemp_all  :', len(templatestr[len(upstreamstr)-10:-len(downstreamstr)]), templatestr[len(upstreamstr)-10:-len(downstreamstr)], 'bandw=', bandw
	if (not currep2==None): print '\ntemp_all  :', len(templatestr[len(upstreamstr)-10:-len(downstreamstr)]), templatestr[len(upstreamstr)-10:-len(downstreamstr)], 'bandw=', bandw
	allsimreads = findTrinucleotideRepeats.myReadTxtFile(curreadfile)

	repregion_len_threhold = len_repPat #3;
	repeatbeforeafter = isupdown
	toleratebeforeafter = 30+isupdown;
	if toleratebeforeafter>len(upstreamstr): toleratebeforeafter = len(upstreamstr)-1
	if toleratebeforeafter>len(downstreamstr): toleratebeforeafter = len(downstreamstr)-1

	if (currep2==None): logging.info("bwamem_w_option=%d, repeatbeforeafter=%d, toleratebeforeafter=%d\n" % (bandwo, repeatbeforeafter, toleratebeforeafter))

	repeatlength = []; replen = len(repregion)
	if seqbp: bp = getAlignment.getBasePair();

	hmmoptions = findTrinucleotideRepeats.getHMMOptions(repeatbeforeafter, repPat, forw_rerv, commonOptions)

	curmaxlen = 0;
	isprint = 0; #1;	
	for j in range(len(allsimreads)-3, 0, -4):
		if mdebug: print "cur line=", j 
		
		if len(allsimreads[j])>3000:
			if curmaxlen<len(allsimreads[j]): curmaxlen=len(allsimreads[j])
			#print 'useUnsymmetricalAlign', j, curmaxlen, len(allsimreads[j]); sys.stdout.flush()

		if ((not commonOptions==None) and len(allsimreads[j])>commonOptions['MaxRep']*len_repPat): 
			logging.error("Error!!!! the sequence is too long!!!"+str(len(allsimreads[j]))+" "+str(j))
			continue;

		querystr = allsimreads[j];
		querystr = string.strip(querystr)

		alignres = UnsymmetricPairAlignment.unsymmetricPairWiseAlignment(templatestr, len(templatestr), querystr, len(querystr), match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandw, isprint).split(';')

		if seqbp:
			querystrR = getAlignment.getComplementary3(bp, querystr)
			alignresR = UnsymmetricPairAlignment.unsymmetricPairWiseAlignment(templatestr, len(templatestr), querystrR, len(querystrR), match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandw, isprint).split(';')

			if int(alignresR[2])>int(alignres[2]):
				querystr = querystrR
				alignres = alignresR
		
		
		if mdebug: print 'done: '
		if mdebug: print 'temp_first', len(templatestr), len(alignres[1]), len(alignres[1][len(upstreamstr)-10:-len(downstreamstr)]), alignres[1][len(upstreamstr)-10:-len(downstreamstr)]
		if mdebug: print 'quey_first', len(querystr), len(alignres[0]), len(alignres[0][len(upstreamstr)-10:-len(downstreamstr)]), alignres[0][len(upstreamstr)-10:-len(downstreamstr)]
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
                                logging.warning("0="+alignres[0])
                                logging.warning("1="+alignres[1])
                                logging.warning("templatestr"+str(len(templatestr))+"="+templatestr)
                                logging.warning("querystr"+str(len(querystr))+"="+querystr)
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
			
	
				if mdebug: print 'temp_rep', len(alignres[1][repeat_start_end[0]:repeat_start_end[1]]), alignres[1][repeat_start_end[0]:repeat_start_end[1]]
				if mdebug: print 'quer_rep', len(alignres[0][repeat_start_end[0]:repeat_start_end[1]]), alignres[0][repeat_start_end[0]:repeat_start_end[1]]
	
				detectregion = alignres[0][repeat_start_end[0]:repeat_start_end[1]].replace('-','');
				if len(detectregion)<repregion_len_threhold and printTwoTail:
					logging.warning("Not enough repeat region is detected")
					logging.warning('temp_rep: %d, %s' % (len(alignres[1][repeat_start_end[0]:repeat_start_end[1]]), alignres[1][repeat_start_end[0]:repeat_start_end[1]]))
					logging.warning('quer_rep: %d, %s' % (len(alignres[0][repeat_start_end[0]:repeat_start_end[1]]), alignres[0][repeat_start_end[0]:repeat_start_end[1]]))
					continue;
	
				if mdebug: print len(detectregion), repeat_start_end, detectregion
				
				if commonOptions==None or (len(detectregion)<commonOptions['MaxRep']*len_repPat):
					if isRemInDel>0:
						newstr, pre0, predstats = findTrinucleotideRepeats.getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, hmmoptions, detectregion, commonOptions)

					else:
						newstr, pre0, predstats = myHMM.hmmpred(detectregion, repPat, forw_rerv, hmmoptions, commonOptions, repeatbeforeafter)
					repeatlength.append(len(newstr)/float(len_repPat)) #3)
				else:
					logging.warning('The sequence is too long: '+str(len(detectregion))+' '+repeatgene+' '+repPat+" "+str(commonOptions['MaxRep'])+" "+str(commonOptions['MaxRep']*len_repPat))

				if repeatlength[-1]>70 and printTwoTail:
					logging.info("More repeat %3d: %s" % (repeatlength[-1], detectregion))
					logging.info("            %3d: %s" % (repeatlength[-1], predstats))
					logging.info("            %3d: %s" % (repeatlength[-1], newstr))
				
	p2, allocr = findTrinucleotideRepeats.get2Peaks(repeatlength, MinSup, commonoptions=commonOptions)
	
	
	return [p2, allocr];

