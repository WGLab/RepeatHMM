
import re;
import os;
import sys;
import string;
import math;
import copy

import numpy;
import peakutils;
import argparse;

import logging

from scipy.stats import norm

import heapq



import getAlignment
import myHMM


from myheader import *

#LOG_FILENAME = 'TrinRepDis.log'
#logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO,filemode='w',format="%(levelname)s: %(message)s")

def myReadTxtFile(filename):
        f = open(filename,'r')

        data = f.readlines();
        while string.strip(data[-1])=="": data = data[:-1];
        f.close();

        return data;

def myWriteTxtFile(mlist, filename):
        f = open(filename, 'w')

        for ml in mlist:
                if ml[-1]=='\n': f.write(ml);
                else: f.write(ml+'\n');

        f.close();

def getDiseseGeneInRefGenomeLocation(hg):
	if hg=='' or hg=='hg38': return getDiseseGeneInRefGenomeLocation_hg38();
	else:
		logging.error('Error: not supported yet')

#from wiki: https://en.wikipedia.org/wiki/Trinucleotide_repeat_disorder
#RefSeq chr:start-end would be used.	
def getDiseseGeneInRefGenomeLocation_hg38():
        gLoc = {}
        #HTT at chr4:3074510-3243960
        gLoc['htt'] = ['chr4','3074510','3243960','3074877','3074939', 'CAG', '+19', '6-35:36-250']                # (CAG)* +1CAACAG 
        #ATN1 at         chr12:   6924463 -  6942321   6936729   6936773
        gLoc['atn1'] = ['chr12', '6924463', '6942321','6936717', '6936773', 'CAG', '+', '6-35:49-88']   # 2(CAGCAA) (CAG)*
        #ar: (Q)* <+ 5non-Q + 7 Q>
        gLoc['ar'] = ['chrX', '67544032', '67730619', '67545318', '67545386', 'CAG', '+22', '9-36:38-62']         # (CAG)* +1CAA
        gLoc['atxn1'] = ['chr6', '16299112', '16761490', '16327636', '16327722', 'CAG', '-14', '6-35:49-88']      # (CAG)* +2(ATGCTG) 11(CTG)
        gLoc['atxn2'] = ['chr12', '111452214', '111599676', '111598951', '111599019', 'CAG', '-13', '14-32:33-77'] # 10Q  (CAG)*
        #atxn3: chr:start-end for RefSeq is chr14:92058552-92106621. chr:start-end for UCSC is chr14:92038652-92106610
        #ATXN3 at chr14:92058552-92106621 
        gLoc['atxn3'] = ['chr14', '92038652', '92106610', '92071011', '92071052', 'CAG', '-14', '12-40:55-86']       # (CAG)* + CTGTTGCTGCTTTTGCTGCTG
        #CACNA1A at chr19:13206443-13506460 
        gLoc['cacna1a'] = ['chr19', '13206443', '13506460', '13207859', '13207897', 'CAG', '-', '4-18:21-30']     #y
        #ATXN7 at chr3:63864557-64003460 
        #could not be found using UCSC genome browser with "simple repeat" on
        gLoc['atxn7'] = ['chr3','63864557', '64003460', '63912686', '63912715', 'CAG', '+7', '7-17:38-120']        # (CAG)* + 3CAG
        #TBP at chr6:           170554333 -  170572870
        gLoc['tbp'] = ['chr6', '170554333', '170572870', '170561899', '170562021', 'CAG', '+19', '25-42:47-63']    # 20*3 (CAG)* +1CAACAG
        #FMR1 at chrX:147911951-147951127 
        gLoc['fmr1'] = ['chrX', '147911951', '147951127', '147912051', '147912110', 'CGG', '+', '6-53:230+/55-200']     # (CGG)* + 1AGG + 9(CGG)
        #aff2: 148500612--148500639: 6 CCG + 3 non-CCG
        #gLoc['aff2'] = ['chrX', '148500619', '149000663', '148500639', '148500692', 'CCG', '+15', '6-35:200+']   # (CCG)* + 1CTG + 2CCG
        gLoc['aff2'] = ['chrX', '148500619', '149000663', '148500606', '148500692', 'CCG', '+15', '6-35:200+']   
        #chr9:69,037,262-69,037,374    FXN at chr9:69035563-69100178
        #could not be found using UCSC genome browser with "simple repeat" on
        gLoc['fxn'] = ['chr9', '69035563', '69100178', '69037287', '69037305', 'GAA', '+5', '7-34:100+']         #xxxy
        #RefSeq for dmpk: chr19:45769709-45782557 
        gLoc['dmpk'] = ['chr19', '45769709', '45782557', '45770205', '45770264', 'CTG', '-20', '5-37:50+']      #y
        #gLoc['sca8'] = ['chr13', '70681345', '70713561', '70139384', '70139428', 'CTG', '15']      #XXXXy
        #gLoc['sca8'] = ['chr13', '70107213', '70139552', '70139351', '70139428', 'CTG', '15']      #Xy  1TTA + 10 CTA + (CTG)*
        #ATXN8OS at       chr13 :  70107213 -  70139753
        gLoc['atxn8os'] = ['chr13', '70107213', '70139753', '70139351', '70139428', 'CTG', '+15', '16-37:110-250']   # 1TTA + 10 CTA + (CTG)*
       
	#https://en.wikipedia.org/wiki/Spinocerebellar_ataxia 
        #PPP2R2B at chr5:           146589505 -  147081520 RefSeq
        gLoc['ppp2r2b'] = ['chr5', '146589505', '147081520', '146878729', '146878759', 'CAG', '-10', '7-28:66-78'] #may contain errors;

        return gLoc;

def get_gLoc(dis_type, gLoc):
        dis_type = dis_type.lower();
        return gLoc[dis_type];

def getValues(mdict):
        valdict = {}
        mdkeys = mdict.keys();
        for mdk in mdkeys:
                if not valdict.has_key(mdict[mdk]):
                        valdict[mdict[mdk]] = 0;
                valdict[mdict[mdk]] += 1;
        valkeys = valdict.keys(); valkeys.sort();
        more = []
        #for i in range(1, len(valkeys)*2/3+1):
        for i in range(1, len(valkeys)+1):
                if valkeys[-i]<4: continue;
                #print -i, valkeys[-i], valdict[valkeys[-i]]
                if valdict[valkeys[-i]]>=2: more.append(valkeys[-i])
        return more;

def findSameValues(mdict, v, mdkeys):
        start = -1; end = -1;
        for mdk_ind in range(len(mdkeys)-1):
                if mdict[mdkeys[mdk_ind]]==mdict[mdkeys[mdk_ind+1]]:
                        if start==-1: start = mdk_ind
                        end = mdk_ind+1
                if (not mdict[mdkeys[mdk_ind]]==mdict[mdkeys[mdk_ind+1]]) and (not start==-1): break;
        return [start, end]

def reviseDict(mdict, v):
        mdkeys = mdict.keys();  mdkeys.sort()
        while True:
                [start, end] = findSameValues(mdict, v, mdkeys)
                if start==-1 or end==-1: break;
                else:
                        i=start; j = end;
                        while True:
                                #print i, j, mdkeys[i], mdkeys[j], mdict[mdkeys[i]], mdict[mdkeys[j]]
                                if i>=j:
                                        print 'Error i is larger than j', i, j, start, end;
                                        break;
                                i += 1;
                                mdict[mdkeys[i]] += 1;
                                if i==j: break;
                                j -= 1;
                                if i==j: break;
                                mdict[mdkeys[j]] += 1;
                                #print i, j, mdkeys[i], mdkeys[j], mdict[mdkeys[i]], mdict[mdkeys[j]]

        return mdict;

def reviseDictAccordingV(mdict):
        more = getValues(mdict)
        for m in more:
                mdict = reviseDict(mdict, m);
        return mdict


def selectFromTwoX(x_index2a, x_index2b, lendict, x_index1):
	if x_index2a==x_index2b: x_index2 = x_index2b
	else:
		neighbor3 = [0,0]; point3 = [x_index2a, x_index2b]
		for pi in range(len(point3)):
			for xi in range(point3[pi]-1, point3[pi]+2):
				if lendict.has_key(xi): neighbor3[pi] += lendict[xi]
		if neighbor3[0]*x_index2a < neighbor3[1]*x_index2b:
			x_index2 = x_index2b
		elif neighbor3[0]*x_index2a > neighbor3[1]*x_index2b:
			x_index2 = x_index2a
		elif x_index2a*2 < x_index2b and neighbor3[0]*0.9 < neighbor3[1]:
			x_index2 = x_index2b
		elif x_index2b*2 < x_index2a and neighbor3[1]*0.9 < neighbor3[0]:
			x_index2 = x_index2a
		else:
			if neighbor3[0]>neighbor3[1]: x_index2 = x_index2a
			elif neighbor3[0]<neighbor3[1]: x_index2 = x_index2b
			else:
				if x_index2a>x_index2b: x_index2 = x_index2a;
				elif x_index2a<x_index2b: x_index2 = x_index2b
				else:
					if math.fabs(x_index1-x_index2a)>math.fabs(x_index1-x_index2b): x_index2 = x_index2a;
					else: x_index2 = x_index2b
	return x_index2

def selectFromTwo(y_index2a, y_index2b, lendict, y_index1, x):
	if y_index2a==y_index2b: y_index2 = y_index2b
	else:
		x2 = selectFromTwoX(x[y_index2a], x[y_index2b], lendict, x[y_index1])
		if x2==x[y_index2a]: y_index2 = y_index2a
		else: y_index2 = y_index2b
	return y_index2
	

def selectOne(iy, curyvalue, lendict, y_index1, indexes, x):
	y_index2a = indexes[iy.index(curyvalue)];
	iy.reverse();
	y_index2b = indexes[len(iy) - iy.index(curyvalue) - 1]

	if y_index2a==y_index1: return y_index2b
	elif y_index2b==y_index1: return y_index2a
	
	return selectFromTwo(y_index2a, y_index2b, lendict, y_index1, x)


def getPeaks2(x, y, lendict, mm, mdebug):
	indexes = peakutils.indexes(numpy.array(y), thres=0.01) - mm
	peak2 = [];
	pnearby = []; 
	
	if mdebug:
		for i in range(len(x)): print ('%d=%d;' % (x[i], y[i])),
		print 'mindexes',
		for i in indexes: print x[i],
		print ''
		print '#x =',
		for i in range(len(x)): print ('%3d' % x[i]),
		print ''
		print 'y  =',
		for i in range(len(x)): print ('%3d' % y[i]),
		print ''

	if len(indexes)>1:
		iy = []
		for i in indexes: iy.append(y[i]);
		ylargest = heapq.nlargest(2, iy);
		y_index1 = indexes[iy.index(ylargest[0])];
		
		p_y_index2 = selectOne(copy.deepcopy(iy), ylargest[1], lendict, y_index1, indexes, x)
		if len(iy)>2:
			ylargest = heapq.nlargest(3, iy);
			p_y_index3 = selectOne(iy, ylargest[2], lendict, y_index1, indexes, x)
			y_index2 = selectFromTwo(p_y_index2, p_y_index3, lendict, y_index1, x)
		else: y_index2 = p_y_index2


		#ylargest = heapq.nlargest(2, iy);
		#y_index1 = indexes[iy.index(ylargest[0])];
		#
		#y_index2a = indexes[iy.index(ylargest[1])];
		#iy.reverse();
		#y_index2b = indexes[len(iy) - iy.index(ylargest[1]) - 1]
		#
		#if y_index2a==y_index2b: y_index2 = y_index2b
		#else:
		#	neighbor3 = [0,0]; point3 = [y_index2a, y_index2b]
		#	for pi in range(len(point3)):
		#		for xi in range(point3[pi]-1, point3[pi]+2):
		#			if lendict.has_key(xi):	neighbor3[pi] += lendict[xi]
		#	if neighbor3[0]>neighbor3[1]: y_index2 = y_index2a
		#	elif neighbor3[0]<neighbor3[1]: y_index2 = y_index2b
		#	else:
		#		if y_index2a>y_index2b: y_index2 = y_index2a;
		#		elif y_index2a<y_index2b: y_index2 = y_index2b
		#		else:
		#			if math.fabs(y_index1-y_index2a)>math.fabs(y_index1-y_index2b): y_index2 = y_index2a;
		#			else: y_index2 = y_index2b

		peak2.append(x[y_index1]); 
		for ni in range(x[y_index1]-2, x[y_index1]+3):
			pnearby.append(ni);
		#if y_index1-1>0: pnearby.append(x[y_index1-1]);
		#if y_index1+1<len(x): pnearby.append(x[y_index1+1]);
		peak2.append(x[y_index2])
		for ni in range(x[y_index2]-2, x[y_index2]+3):
			pnearby.append(ni);
		#if y_index2-1>0: pnearby.append(x[y_index2-1]); 
		#if y_index2+1<len(x): pnearby.append(x[y_index2+1]);
        else:
		if len(indexes)==1:
			logging.debug("Only find one peak");
			peak2.append(x[indexes[0]])
			for ni in range(x[indexes[0]]-2, x[indexes[0]]+3):
				pnearby.append(ni);
			#pnearby.append(x[indexes[0]]);
			#if indexes[0]-1>0: pnearby.append(x[indexes[0]-1]); 
			#if indexes[0]+1<len(x): pnearby.append(x[indexes[0]+1]);
		else:
			logging.debug("Cannot find a peak");

	if len(peak2)>1 and y[peak2.index(peak2[1])] > [y[peak2.index(peak2[0])]]:
		peak2[0], peak2[1] = peak2[1], peak2[0]
	
	p2len = len(peak2)
	twotails = [];
	if (x[0] not in pnearby): twotails.append(x[0]);
	if (x[len(x)-1] not in pnearby): twotails.append(x[len(x)-1]);
        #if (x[0] not in peak2) and (x[1] not in peak2):
        #       if lendict[x[0]]==lendict[x[1]]: twotails.append(x[0]);
        #       else: twotails.append(x[0]);twotails.append(x[1]);
        #if (x[len(x)-2] not in peak2) and (x[len(x)-1] not in peak2):
        #       if lendict[x[len(x)-2]]==lendict[x[len(x)-1]]:  twotails.append(x[len(x)-1]);
        #       else: twotails.append(x[len(x)-2]); twotails.append(x[len(x)-1]);
	
	for ti in twotails:
		for i in range(p2len):
			if selectFromTwoX(ti, peak2[i], lendict, -1)==ti:
			#if lendict[ti] > lendict[peak2[i]]: 
				peak2.insert(i, ti); break;

	if len(peak2)>2: peak2 = peak2[:2];

	return peak2	

def get2Peaks(lengd):
	#prnucleotideRepeats.pyint lengd

	pvalue = 0.05; ratio = 0.4;

	mdebug = True; mdebug = False;
	
	lendict = {};
	for l in lengd:
		l = int(l+0.5)
		if not lendict.has_key(l): lendict[l] = 0;
		lendict[l] += 1;

	ldkeys = lendict.keys(); ldkeys.sort();
	
	allocr = '';
	for ldk in ldkeys:
		allocr += ('%d:%d, ' % (ldk, lendict[ldk]))
	logging.info(allocr)

	if len(ldkeys)<1: return [[0], allocr]
	elif len(ldkeys)<2: return [[ldkeys[0]], allocr]


	len3dict = {}
	for i in range(len(ldkeys)):
                curlk = ldkeys[i]
                if i==0:
                        if curlk==0: len3dict[curlk] = 0
                        else:
                                len3dict[curlk] = lendict[curlk]
                                if lendict.has_key(curlk+1): len3dict[curlk] += lendict[curlk+1]
                                else: len3dict[curlk+1] = 0
                elif i==len(ldkeys)-1:
                        if curlk-1==0: len3dict[curlk] = lendict[curlk]
                        else:
                                len3dict[curlk] = lendict[curlk]
                                if lendict.has_key(curlk-1): len3dict[curlk] += lendict[curlk-1]
                                else: len3dict[curlk-1] = 0
                else:
                        if curlk-1==0:
                                len3dict[curlk] = lendict[curlk]
                                if lendict.has_key(curlk+1): len3dict[curlk] += lendict[curlk+1]
                                else: len3dict[curlk+1] = 0
                        else:
                                len3dict[curlk] = lendict[curlk]
                                if lendict.has_key(curlk-1): len3dict[curlk] += lendict[curlk-1]
                                else: len3dict[curlk-1] = 0
                                if lendict.has_key(curlk+1): len3dict[curlk] += lendict[curlk+1]
                                else: len3dict[curlk+1] = 0
	'''
	for i in range(len(ldkeys)):
                if i==0:
	                if ldkeys[i]==0: len3dict[ldkeys[i]] = 0; #lendict[ldkeys[i]]
                        else: len3dict[ldkeys[i]] = lendict[ldkeys[i]] + lendict[ldkeys[i+1]]
                elif i==len(ldkeys)-1:
                        if ldkeys[i-1]==0: len3dict[ldkeys[i]] = lendict[ldkeys[i]]
                        else: len3dict[ldkeys[i]] = lendict[ldkeys[i]] + lendict[ldkeys[i-1]]
                else:
                        if ldkeys[i-1]==0: len3dict[ldkeys[i]] = lendict[ldkeys[i]] + lendict[ldkeys[i+1]]
                        else: len3dict[ldkeys[i]] = lendict[ldkeys[i-1]] + lendict[ldkeys[i]] + lendict[ldkeys[i+1]]
	'''
	
	len3dict = reviseDictAccordingV(len3dict)
	
	x = []; yo = [];
	for ldk in ldkeys:
		x.append(ldk); 
		yo.append(len3dict[ldk]);

	peak2 = getPeaks2(x, yo, lendict, 0, mdebug)	


	'''
	lendict = {};
	for l in lengd:
		l = int(l+0.5)
		if not lendict.has_key(l): lendict[l] = 0;
		lendict[l] += 1;

	lendict = reviseDictAccordingV(lendict)
	#print lendict

	x = []; yo = []; allocr = ''; 
	yo = [0]; 
	ldkeys = lendict.keys(); ldkeys.sort();
	for ldk in ldkeys:
		x.append(ldk); yo.append(lendict[ldk]);
		allocr += ('%d:%d, ' % (ldk, lendict[ldk]))

	#print lendict	
	logging.info(allocr)
	
	if len(ldkeys)<1: return [[0], allocr]
	elif len(ldkeys)<2: return [[yo[1]], allocr]

	yo.append(0); 
	
	pea1k2 = getPeaks2(x, yo[1:-1], lendict, 0, mdebug)	
	
	y2a = []
	for i in range(len(yo)-1): y2a.append(yo[i]+yo[i+1])
	pea2k2 = getPeaks2(x, y2a, lendict, 1, mdebug)

	pea1knearby = []; 
	for p1 in pea1k2:
		for ni in range(p1-2, p1+3):
			pea1knearby.append(ni);
		#pea1knearby.append(p1);
		#p1ind = x.index(p1);
		#if p1ind-1>0: pea1knearby.append(x[p1ind-1]);
		#if p1ind+1<len(x): pea1knearby.append(x[p1ind+1]);

	pea2k2new = []
	for p2 in pea2k2:
		p2ind = x.index(p2);
		nb3 = [p2]
		if p2ind-1>0: nb3.append(x[p2ind-1]);
		if p2ind+1<len(x): nb3.append(x[p2ind+1]);
		curmx = nb3[-1]  #p2;
		for curi in nb3: 
			if lendict[curi]>lendict[curmx]: curmx = curi
		pea2k2new.append(curmx)

	notinpea1k = []
	for p2 in pea2k2new:
		if p2 not in pea1knearby: notinpea1k.append(p2);

	#print pea1k2, notinpea1k, pea2k2new 

	#if len(pea1k2)==0: pea1k2 = pea2k2new
	#else:
	for p2 in notinpea1k:
			has_insert = False;
			for i in range(len(pea1k2)):
				x_index_find = selectFromTwoX(p2, pea1k2[i], lendict, -1)
				if x_index_find == p2:
					pea1k2.insert(i, p2); has_insert = True; break;
				#if lendict[p2] > lendict[pea1k2[i]]:
				#	 pea1k2.insert(i, p2); has_insert = True; break;
			if not has_insert: pea1k2.append(p2)

	peak2 = pea1k2
	'''
	if len(peak2)>2: peak2 = peak2[:2]

	if mdebug:
		for i in range(len(x)): print ('%d=%d;' % (x[i], lendict[x[i]])),
		print 'mPeak', peak2
		print '#xF=',
		for i in range(len(x)): print ('%3d' % x[i]),
		print ''
		print 'yF =',
		for i in range(len(x)): print ('%3d' % lendict[x[i]]),
		print ''

	ally = [lendict[ldkeys[0]], lendict[ldkeys[-1]]]
	if len(ldkeys)>1: last2 = -2;
	else: last2 = -1;
	for ik in range(ldkeys[1], ldkeys[last2]+1):
		if (ik in peak2) and len(ldkeys)>1: continue;
		if lendict.has_key(ik): ally.append(lendict[ik]);
		else: ally.append(0)

	#mmean = numpy.mean(numpy.array(yo[1:-1]));
	#mstd = numpy.std(numpy.array(yo[1:-1]));
	mmean = numpy.mean(numpy.array(ally));
	mstd = numpy.std(numpy.array(ally));
	if mdebug:
		print ('<%.2f,%.2f/%d>' % (mmean, mstd, len(yo[1:-1]))),
		for pf in peak2:
			print (' %d=%.6f:%.3f; ' % (pf, 1-norm.cdf((lendict[pf]-mmean)/mstd), (lendict[pf]-mmean)/mstd)), 

	if len(peak2)==2:  # pvalue = 0.05; ratio = 0.55
		mstr = 'In Peak:'
		for pf in peak2:
			mstr += ('for %d, pvalue=%.6f; ' % (pf, 1-norm.cdf((lendict[pf]-mmean)/mstd)))
		if not lendict.has_key(peak2[0]): mstr += (' No %d' % peak2[0]);
		elif not lendict.has_key(peak2[1]): mstr += (' No %d' % peak2[1]);
		else:  mstr += (' ratio=%.3f(%d/%d)' % (lendict[peak2[1]]/float(lendict[peak2[0]]), lendict[peak2[1]], lendict[peak2[0]]))
		#logging.info(mstr)

		if mdebug: 
			print ('%.3f=%d/%d: ' % (lendict[peak2[1]]/float(lendict[peak2[0]]), lendict[peak2[1]], lendict[peak2[0]]) ), 
		#if lendict[peak2[1]]/float(lendict[peak2[0]])<ratio or 1-norm.cdf((lendict[peak2[1]]-mmean)/mstd)>pvalue:
		#if 1-norm.cdf((lendict[peak2[1]]-mmean)/mstd)>pvalue:
		sum3_0 = lendict[peak2[0]]; sum3_1 = lendict[peak2[1]];
		if lendict.has_key(peak2[0]-1) and lendict.has_key(peak2[1]-1):
			sum3_0 += lendict[peak2[0]-1]
			sum3_1 += lendict[peak2[1]-1]
		if lendict.has_key(peak2[0]+1) and lendict.has_key(peak2[1]+1):
			sum3_0 += lendict[peak2[0]+1]
			sum3_1 += lendict[peak2[1]+1]
		
		mstr += (' >>> %.3f(%d/%d)' % (sum3_1*peak2[1]/float(sum3_0*peak2[0]), sum3_1*peak2[1], sum3_0*peak2[0]))
		logging.info(mstr)
		
		if sum3_1*peak2[1]/float(sum3_0*peak2[0])<0.1: 
			peak2 = peak2[:1]

	for curp_ind in range(len(peak2)):
                curp = peak2[curp_ind]
                pin = ldkeys.index(curp)
                nb3 = [curp]
                if pin-1>0: nb3.append(ldkeys[pin-1])
                if pin+1<len(ldkeys): nb3.append(ldkeys[pin+1])
                curmx = nb3[-1]
                for curi in nb3:
                        if lendict[curi]>lendict[curmx]: curmx = curi
                peak2[curp_ind] = curmx

	peak2.sort();
	if len(peak2)>1:
                if len3dict[peak2[0]]<len3dict[peak2[1]]*0.6:
                        peak2 = [peak2[1], peak2[1]]


	if mdebug: print 'mPeak', peak2

	return [peak2, allocr[:-1]];

#def getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, queryrep, match=10, mismatch=-9, gap_in_perf=-2, gap_in_read=-13, gap_before_after = -1):
def getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, queryrep):
	unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, queryrep[repeatbeforeafter:(len(queryrep)-repeatbeforeafter)], forw_rerv, match, mismatch, gap_in_perf, gap_in_read, gap_before_after);
	if repeatbeforeafter>0: unewstr =  queryrep[:repeatbeforeafter]+ unewstr + queryrep[(len(queryrep)-repeatbeforeafter):]
	newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, repeatbeforeafter)
	
	return [newstr, pre0, predstats]
	

def getGene(repeatgene, chr, gene_start_end, unique_file_id, analysis_file_id):
        alignfolder = 'align/'
	if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

        #unique_file_id = simulation_file_id + analysis_file_id

        fastafile = alignfolder + repeatgene + unique_file_id +'.fasta'
        #get_alg_cmd = 'samtools faidx ./hg38/Homo_sapiens_assembly38.fasta '+ chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+fastafile

        get_alg_cmd = 'samtools faidx '+hg38_reference_and_index+'/hg38.fa '+ chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+fastafile
        print get_alg_cmd

        os.system(get_alg_cmd)
        if not os.path.isfile(fastafile):
                logging.error('Cannot produce '+fastafile+' for '+repeatgene)
                sys.exit(1)
        fastadata = myReadTxtFile(fastafile)
        mfadata = ''
        for li in fastadata:
            if li[0]=='>': continue;
            if li[-1]=='\n': li = li[:-1]
            mfadata = mfadata + li

	os.system('rm '+fastafile)

	return mfadata.upper()

def getRepeatForGivenGene2(chr, repeatgene, gene_start_end, repeat_orig_start_end, bamfile, repPat, forw_rerv, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id):
        #print repeatgene,
        alignfolder = 'align/'
        if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

        ref_repeat = (repeat_orig_start_end[1]-repeat_orig_start_end[0]+1)/3.0
        alignfile = alignfolder + repeatgene + unique_file_id +'.alignment.txt'
        get_alg_cmd = 'samtools view '+bamfile+' ' + chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+alignfile
        logging.info('Running '+get_alg_cmd)
        os.system(get_alg_cmd);
        logging.info('Produced ' + alignfile + ' done!');

        if not os.path.isfile(alignfile):
                logging.error('Cannot produce '+alignfile+' for '+repeatgene)
                sys.exit(1)
        aligndata = myReadTxtFile(alignfile)
        os.system('rm '+alignfile)

        mfadata = getGene(repeatgene, chr, gene_start_end, unique_file_id, analysis_file_id)

        covermorebeforeafter = 30
        repregion_len_threhold = 3;
        repeatbeforeafter = isupdown - isExtend
        repeat_start_end = [repeat_orig_start_end[0], repeat_orig_start_end[1]]
        #repeat_start_end[0] -= isExtend; repeat_start_end[1] += isExtend;
        #if isExtend>0 and repeat_start_end[0]<1: repeat_start_end[0]=1
	
        simplebeforeafter = 3;

        check = False; #True;
        wrongalign = 0;
        if check: repeat_beforeafter = [];

        repeats = [];
        for line in aligndata:
                if check:
                   beforematch = {}; aftermatch = {}; inmatch = 0; allmatch = 0; neighbeforeafter = 200;
                   beforematch[50] = 0; aftermatch[50] = 0
                   beforematch[100] = 0; aftermatch[100] = 0
                   beforematch[150] = 0; aftermatch[150] = 0
                   beforematch[200] = 0; aftermatch[200] = 0
                   beforematch[250] = 0; aftermatch[250] = 0
                   beforematch[300] = 0; aftermatch[300] = 0

                lsp = line.split('\t')
                cchr = lsp[2]
                pos = int(lsp[3])
                aligninfo = lsp[5]
                aainfo = lsp[9]
                #qualifyinfo = lsp[12]
                #print 'qualifyinfo', qualifyinfo[:100]

                if pos > repeat_start_end[0] - covermorebeforeafter:
                        wrongalign += 1;
                        continue;
                        #logging.error('The start pos in ref Genome is greater than the start position of repeats' + str(pos) +' ' + str(repeat_start_end[0]));
                if not cchr==chr:  logging.error('Not same ' + cchr +' ' + chr); continue;

                numreg = re.compile('\d+')
                numinfo = numreg.findall(aligninfo)

                mdireg = re.compile('[MIDNSHPX=]{1}')
                mdiinfo = mdireg.findall(aligninfo)

                if not len(numinfo)==len(mdiinfo):
                        logging.error('Num is equal to mid' +str(len(numinfo)) + ' '+ str(len(mdiinfo))); continue;

                queryind = 0;
                queryrep = '';
                longer = False;
                if check:
                   totalbefore = 0; totalafter = 0;  matchinfo = '';
                   neighmatch = 0; neighref=''; neightest = ''

                for n1ind in range(len(numinfo)):
                        n1 = int(numinfo[n1ind])
                        mdi = mdiinfo[n1ind];

                        for n1i in range(n1):
                                if check:
                                   if totalbefore<repeat_start_end[0]-pos:
                                      totalbefore=repeat_start_end[0]-pos;
                                   if totalafter<pos-repeat_start_end[1]:
                                      totalafter=pos-repeat_start_end[1]

                                qrepadd = False;
                                if mdi=='M':
                                        if check:
                                           faind = pos - (repeat_start_end[0] -  neighbeforeafter)
                                           if (faind>=0 and pos-repeat_start_end[0]<0) or \
                                              (pos-repeat_start_end[1]>=0 and pos-(repeat_start_end[1]+neighbeforeafter)<=0):
                                                if mfadata[faind] == aainfo[queryind]: neighmatch += 1;
                                                neighref = neighref + mfadata[faind]
                                                neightest = neightest + aainfo[queryind]

                                        pos = pos + 1;
                                        queryind = queryind + 1;
                                        qrepadd = True;

                                        if check:
                                           allmatch += 1;
                                           if pos-1 < repeat_start_end[0]:
                                                bef = repeat_start_end[0] - pos + 1
                                                bmkeys = beforematch.keys(); bmkeys.sort();
                                                for bmk in bmkeys:
                                                    if bmk>=bef: beforematch[bmk] += 1;
                                           elif pos-1 > repeat_start_end[1]:
                                                aft = pos-1 - repeat_start_end[1]
                                                afkeys = aftermatch.keys(); afkeys.sort();
                                                for afk in afkeys:
                                                    if afk >= aft: aftermatch[afk] += 1;
                                           else: inmatch += 1;

                                elif mdi =='I':
                                        qrepadd = True;
                                        queryind = queryind + 1;
                                elif mdi == 'D':
                                        pos = pos + 1;
                                elif mdi == 'S':
                                        queryind = queryind + 1;
                                        qrepadd = True;
                                elif mdi == 'H':
                                        pass;
                                elif mdi == 'P':
                                        pass;
                                else:
                                        logging.warning('Warning unknow CIGAR element ' + str(n1) + ' ' + mdi)
                                if qrepadd:
                                        #if pos-1 >= repeat_start_end[0]-repeatbeforeafter and pos-1 <= repeat_start_end[1]+repeatbeforeafter:
                                        if pos-1 >= repeat_start_end[0]-simplebeforeafter and pos-1 <= repeat_start_end[1]+simplebeforeafter:
                                                queryrep = queryrep + aainfo[queryind-1]
                                if check and pos-1 >= repeat_start_end[0]-repeatbeforeafter and pos-1 <= repeat_start_end[1]+repeatbeforeafter: matchinfo += mdi
                        #if pos-1 > repeat_start_end[1] + covermorebeforeafter: longer = True;
                        if pos-1 > repeat_start_end[1]+simplebeforeafter: longer = True;

                if check and len(queryrep)>=repregion_len_threhold:
                   bmkeys = beforematch.keys(); bmkeys.sort(); befstr = '';
                   for bmk in bmkeys:
                       befstr += (str(bmk)+'='+str(beforematch[bmk])+';');
                   afkeys = aftermatch.keys(); afkeys.sort(); aftstr = '';
                   for afk in afkeys:
                       aftstr += (str(afk)+'='+str(aftermatch[afk])+';');

                   print repeatgene, repeat_start_end, len(queryrep), queryrep, gene_start_end

                   newstr = '';
                   print repPat, ':', 'MI', matchinfo
                   if len(queryrep)<1000:
                        newstr, pre0, predstats = getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, queryrep)

                        #unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, queryrep[repeatbeforeafter:(len(queryrep)-repeatbeforeafter)], forw_rerv);
                        #if repeatbeforeafter>0: unewstr =  queryrep[:repeatbeforeafter]+ unewstr + queryrep[(len(queryrep)-repeatbeforeafter):]
                        #newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, repeatbeforeafter)

                        #newstr, pre0, predstats = myHMM.hmmpred(queryrep, repPat, forw_rerv, repeatbeforeafter)
                   print ('Gene %s(%s:%d-%d) Match for repeat(%s%s):INmatch=%d/%d(%d) test_rep=%d; beforematch:%s(%d) aftermatch:%s(%d)' % (repeatgene, chr, repeat_start_end[0], repeat_start_end[1], repPat, forw_rerv[0], inmatch, len(queryrep), neighmatch, len(newstr)/3, befstr, totalbefore, aftstr, totalafter)), longer
                   #print neighref;
                   #print neightest;

                   if longer: repeat_beforeafter.append([len(newstr)/3-simplebeforeafter*2/3, totalbefore, totalafter, neighmatch])

                if len(queryrep)>=repregion_len_threhold: repeats.append([longer, queryrep, lsp[0]])

        if check:
           for rbf in repeat_beforeafter:
                print ('%6d %6d %6d %6d' % (rbf[0], rbf[1], rbf[2], rbf[3]))

        rptrue = []; rpfalse = []; orignial = [];
        for currep in repeats:
                #print chr, repeatgene, repPat, len(currep[1]),
                newstr = currep[1]

                pre0 = 0; predstats=''
                if len(newstr)<1000:
                        #newstr = getAlignment.correctSeq(repPat, currep[1], forw_rerv);
                        if isAlign>0:
                                pass #newstr, pre0, predstats = getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, currep[1])

                                #unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, currep[1][repeatbeforeafter:(len(currep[1])-repeatbeforeafter)], forw_rerv);
                                #if repeatbeforeafter>0: unewstr =  currep[1][:repeatbeforeafter]+ unewstr + currep[1][(len(currep[1])-repeatbeforeafter):]
                                #
                                ##newstr, pre0, predstats = myHMM.hmmpred(currep[1], repPat, forw_rerv, repeatbeforeafter)
                                #newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, repeatbeforeafter)


                        else:
                                pass #newstr, pre0, predstats = myHMM.hmmpred(newstr, repPat, forw_rerv, repeatbeforeafter)
                else: logging.warning('The sequence is too long: '+str(len(newstr))+' '+chr+' '+repeatgene+' '+repPat+' '+str(currep[0])+' reads name:'+currep[2])
                orignial.append([currep[1], pre0, predstats]);
                currep[1] = newstr
                #if len(currep[1])==0: continue;
                if currep[0]: #repPat, forw_rerv
                        #rptrue.append(len(newstr)/3.0-repeatbeforeafter/3);
                        #if isHMM: rptrue.append(len(currep[1])/3.0);
                        rptrue.append(len(currep[1])/3.0-simplebeforeafter*2/3);
                else:
                        #rpfalse.append(len(newstr)/3.0-repeatbeforeafter/3);
                        #if isHMM: rpfalse.append(len(currep[1])/3.0);
                        rpfalse.append(len(currep[1])/3.0-simplebeforeafter*2/3);

        rptrue.sort(); rpfalse.sort()
        trstr = 'true ' + str(len(rptrue)) + ' [';
        for rpt in rptrue:
                trstr = trstr + ('%.0f,' % rpt)
        trstr = trstr[:-1] + ']'
        logging.debug(trstr)

        #print repeatgene, repPat, rptrue
        p2, allocr = get2Peaks(rptrue)

        if len(rpfalse)>0:
                flstr = 'fals ' + str(len(rpfalse)) + ' ['
                for rpf in rpfalse:
                        flstr = flstr + ('%.0f,' % rpf)
                flstr = flstr[:-1] + ']'
                logging.debug(flstr);

        logging.info('ref_repeat ' + ('%.0f' % ref_repeat) +'\t'+repPat+'\t'+forw_rerv);

        '''
        for currep_ind in range(len(repeats)):
                currep = repeats[currep_ind]

                aaprinindex = -1;
                if not (currep[0]): aaprinindex = 300

                logging.debug('\t'+str(currep[0]) + ' o:' + str(len(orignial[currep_ind][0]))  +'\t'+ orignial[currep_ind][0][:aaprinindex]);
                prestr = '';
                for i in range(orignial[currep_ind][1]): prestr += ' ';
                #logging.debug('\t'+str(currep[0]) + ':' + str(len(orignial[currep_ind][0]))  +'\t'+prestr+ orignial[currep_ind][2]);
                #logging.debug('\t'+str(currep[0]) + ':' + str(len(orignial[currep_ind][0]))  +'\t'+orignial[currep_ind][2]);
                logging.debug('\t'+str(currep[0]) + ' p:' + str(len(currep[1])) +'\t' + prestr+ (currep[1][:aaprinindex]))
        '''
        #p2, allocr = get2Peaks(rptrue)

        return [repeatgene, ref_repeat, p2, allocr, len(rptrue), len(rpfalse)+wrongalign]


def getRepeatForGivenGene(chr, repeatgene, gene_start_end, repeat_orig_start_end, bamfile, repPat, forw_rerv, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id):
        #print repeatgene,
        alignfolder = 'align/'	
	if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

        ref_repeat = (repeat_orig_start_end[1]-repeat_orig_start_end[0]+1)/3.0
        #repeat_start_end[1] += 1;

	##unique_file_id = simulation_file_id + analysis_file_id
        alignfile = alignfolder + repeatgene + unique_file_id +'.alignment.txt'
        get_alg_cmd = 'samtools view '+bamfile+' ' + chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+alignfile
        logging.info('Running '+get_alg_cmd)
        os.system(get_alg_cmd);
        logging.info('Produced ' + alignfile + ' done!');

        if not os.path.isfile(alignfile):
                logging.error('Cannot produce '+alignfile+' for '+repeatgene)
                sys.exit(1)
        aligndata = myReadTxtFile(alignfile)
	os.system('rm '+alignfile)

	mfadata = getGene(repeatgene, chr, gene_start_end, unique_file_id, analysis_file_id)
	
        #fastafile = alignfolder + repeatgene+'.1.fasta'
        #get_alg_cmd = 'samtools faidx ./hg38/Homo_sapiens_assembly38.fasta '+ chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+fastafile
        #os.system(get_alg_cmd)
        #if not os.path.isfile(fastafile):
        #        logging.error('Cannot produce '+fastafile+' for '+repeatgene)
        #        sys.exit(1)
        #fastadata = myReadTxtFile(fastafile)
        #mfadata = ''
	#for li in fastadata:
        #    if li[0]=='>': continue;
        #    if li[-1]=='\n': li = li[:-1]
        #    mfadata = mfadata + li

        covermorebeforeafter = 30
        repregion_len_threhold = 3;
        repeatbeforeafter = isupdown - isExtend
        repeat_start_end = [repeat_orig_start_end[0], repeat_orig_start_end[1]]
        repeat_start_end[0] -= isExtend; repeat_start_end[1] += isExtend;
        if isExtend>0 and repeat_start_end[0]<1: repeat_start_end[0]=1

        check = False; #True;
	wrongalign = 0;
        if check: repeat_beforeafter = []; 

        repeats = [];
        for line in aligndata:
                if check:
                   beforematch = {}; aftermatch = {}; inmatch = 0; allmatch = 0; neighbeforeafter = 200;
                   beforematch[50] = 0; aftermatch[50] = 0
                   beforematch[100] = 0; aftermatch[100] = 0
                   beforematch[150] = 0; aftermatch[150] = 0
                   beforematch[200] = 0; aftermatch[200] = 0
                   beforematch[250] = 0; aftermatch[250] = 0
                   beforematch[300] = 0; aftermatch[300] = 0

                lsp = line.split('\t')
                cchr = lsp[2]
                pos = int(lsp[3])
                aligninfo = lsp[5]
                aainfo = lsp[9]
                #qualifyinfo = lsp[12]
                #print 'qualifyinfo', qualifyinfo[:100]

                if pos > repeat_start_end[0] - covermorebeforeafter: 
			wrongalign += 1;
                        continue;  
                        #logging.error('The start pos in ref Genome is greater than the start position of repeats' + str(pos) +' ' + str(repeat_start_end[0]));
                if not cchr==chr:  logging.error('Not same ' + cchr +' ' + chr); continue;

                numreg = re.compile('\d+')
                numinfo = numreg.findall(aligninfo)

                mdireg = re.compile('[MIDNSHPX=]{1}')
                mdiinfo = mdireg.findall(aligninfo)

                if not len(numinfo)==len(mdiinfo):
                        logging.error('Num is equal to mid' +str(len(numinfo)) + ' '+ str(len(mdiinfo))); continue;

                queryind = 0;
                queryrep = '';
                longer = False;
                if check: 
                   totalbefore = 0; totalafter = 0;  matchinfo = ''; 
                   neighmatch = 0; neighref=''; neightest = ''

                for n1ind in range(len(numinfo)):
                        n1 = int(numinfo[n1ind])
                        mdi = mdiinfo[n1ind];

                        for n1i in range(n1):
                                if check:
                                   if totalbefore<repeat_start_end[0]-pos: 
                                      totalbefore=repeat_start_end[0]-pos;
                                   if totalafter<pos-repeat_start_end[1]: 
                                      totalafter=pos-repeat_start_end[1]

                                qrepadd = False;
                                if mdi=='M':
                                        if check:
                                           faind = pos - (repeat_start_end[0] -  neighbeforeafter)
                                           if (faind>=0 and pos-repeat_start_end[0]<0) or \
                                              (pos-repeat_start_end[1]>=0 and pos-(repeat_start_end[1]+neighbeforeafter)<=0):
                                                if mfadata[faind] == aainfo[queryind]: neighmatch += 1;
                                                neighref = neighref + mfadata[faind]
                                                neightest = neightest + aainfo[queryind]

                                        pos = pos + 1;
                                        queryind = queryind + 1;
                                        qrepadd = True;
                                 
                                        if check: 
                                           allmatch += 1;
                                           if pos-1 < repeat_start_end[0]:
                                                bef = repeat_start_end[0] - pos + 1
                                                bmkeys = beforematch.keys(); bmkeys.sort();
                                                for bmk in bmkeys:
                                                    if bmk>=bef: beforematch[bmk] += 1;
                                           elif pos-1 > repeat_start_end[1]:
                                                aft = pos-1 - repeat_start_end[1]
                                                afkeys = aftermatch.keys(); afkeys.sort();
                                                for afk in afkeys:
                                                    if afk >= aft: aftermatch[afk] += 1;
                                           else: inmatch += 1;            

                                elif mdi =='I':
                                        qrepadd = True;
                                        queryind = queryind + 1;
                                elif mdi == 'D':
                                        pos = pos + 1;
                                elif mdi == 'S':
                                        queryind = queryind + 1;
                                        qrepadd = True;
                                elif mdi == 'H':
                                        pass;
                                elif mdi == 'P':
                                        pass;
                                else:
                                        logging.warning('Warning unknow CIGAR element ' + str(n1) + ' ' + mdi)
                                if qrepadd:
                                        if pos-1 >= repeat_start_end[0]-repeatbeforeafter and pos-1 <= repeat_start_end[1]+repeatbeforeafter:
                                                queryrep = queryrep + aainfo[queryind-1]
                                if check and pos-1 >= repeat_start_end[0]-repeatbeforeafter and pos-1 <= repeat_start_end[1]+repeatbeforeafter: matchinfo += mdi
                        if pos-1 > repeat_start_end[1] + covermorebeforeafter: longer = True;

                if check and len(queryrep)>=repregion_len_threhold:
                   bmkeys = beforematch.keys(); bmkeys.sort(); befstr = '';
                   for bmk in bmkeys:
                       befstr += (str(bmk)+'='+str(beforematch[bmk])+';');           
                   afkeys = aftermatch.keys(); afkeys.sort(); aftstr = '';
                   for afk in afkeys:
                       aftstr += (str(afk)+'='+str(aftermatch[afk])+';');
	
                   print repeatgene, repeat_start_end, len(queryrep), queryrep, gene_start_end

                   newstr = ''; 
                   print repPat, ':', 'MI', matchinfo 
                   if len(queryrep)<1000:
                        newstr, pre0, predstats = getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, queryrep)

                        #unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, queryrep[repeatbeforeafter:(len(queryrep)-repeatbeforeafter)], forw_rerv);
                        #if repeatbeforeafter>0: unewstr =  queryrep[:repeatbeforeafter]+ unewstr + queryrep[(len(queryrep)-repeatbeforeafter):]
                        #newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, repeatbeforeafter)

                        #newstr, pre0, predstats = myHMM.hmmpred(queryrep, repPat, forw_rerv, repeatbeforeafter)
                   print ('Gene %s(%s:%d-%d) Match for repeat(%s%s):INmatch=%d/%d(%d) test_rep=%d; beforematch:%s(%d) aftermatch:%s(%d)' % (repeatgene, chr, repeat_start_end[0], repeat_start_end[1], repPat, forw_rerv[0], inmatch, len(queryrep), neighmatch, len(newstr)/3, befstr, totalbefore, aftstr, totalafter)), longer
                   #print neighref; 
                   #print neightest;

                   if longer: repeat_beforeafter.append([len(newstr)/3, totalbefore, totalafter, neighmatch])

                if len(queryrep)>=repregion_len_threhold: repeats.append([longer, queryrep, lsp[0]])

        if check:
           for rbf in repeat_beforeafter:
                print ('%6d %6d %6d %6d' % (rbf[0], rbf[1], rbf[2], rbf[3]))

        rptrue = []; rpfalse = []; orignial = [];  
        for currep in repeats:
                #print chr, repeatgene, repPat, len(currep[1]), 
                newstr = currep[1]
                
                pre0 = 0; predstats=''
                if len(newstr)<1000: 
			#newstr = getAlignment.correctSeq(repPat, currep[1], forw_rerv);
			if isAlign>0:
				newstr, pre0, predstats = getUnsymAlignAndHMM(repPat, forw_rerv, repeatbeforeafter, currep[1])

				#unewstr = getAlignment.myUnsymmetricPairAlignment(repPat, currep[1][repeatbeforeafter:(len(currep[1])-repeatbeforeafter)], forw_rerv);
				#if repeatbeforeafter>0: unewstr =  currep[1][:repeatbeforeafter]+ unewstr + currep[1][(len(currep[1])-repeatbeforeafter):]
				#
				##newstr, pre0, predstats = myHMM.hmmpred(currep[1], repPat, forw_rerv, repeatbeforeafter)
				#newstr, pre0, predstats = myHMM.hmmpred(unewstr, repPat, forw_rerv, repeatbeforeafter)
		

			else:
				newstr, pre0, predstats = myHMM.hmmpred(newstr, repPat, forw_rerv, repeatbeforeafter)
                else: logging.warning('The sequence is too long: '+str(len(newstr))+' '+chr+' '+repeatgene+' '+repPat+' '+str(currep[0])+' reads name:'+currep[2])
                orignial.append([currep[1], pre0, predstats]);
                currep[1] = newstr 
                #if len(currep[1])==0: continue;
                if currep[0]: #repPat, forw_rerv 
                        #rptrue.append(len(newstr)/3.0-repeatbeforeafter/3);
                        #if isHMM: rptrue.append(len(currep[1])/3.0);
                        rptrue.append(len(currep[1])/3.0);
                else: 
                        #rpfalse.append(len(newstr)/3.0-repeatbeforeafter/3);
                        #if isHMM: rpfalse.append(len(currep[1])/3.0);
                        rpfalse.append(len(currep[1])/3.0);

        rptrue.sort(); rpfalse.sort()
        trstr = 'true ' + str(len(rptrue)) + ' [';
        for rpt in rptrue:
                trstr = trstr + ('%.0f,' % rpt)
        trstr = trstr[:-1] + ']'
	logging.debug(trstr)

	#print repeatgene, repPat, rptrue
	p2, allocr = get2Peaks(rptrue)

	if len(rpfalse)>0:
                flstr = 'fals ' + str(len(rpfalse)) + ' ['
                for rpf in rpfalse:
                        flstr = flstr + ('%.0f,' % rpf)
                flstr = flstr[:-1] + ']'
                logging.debug(flstr);

        logging.info('ref_repeat ' + ('%.0f' % ref_repeat) +'\t'+repPat+'\t'+forw_rerv);

        for currep_ind in range(len(repeats)):
                currep = repeats[currep_ind]
             
                aaprinindex = -1;
                if not (currep[0]): aaprinindex = 300

                logging.debug('\t'+str(currep[0]) + ' o:' + str(len(orignial[currep_ind][0]))  +'\t'+ orignial[currep_ind][0][:aaprinindex]);
                prestr = ''; 
                for i in range(orignial[currep_ind][1]): prestr += ' ';
                #logging.debug('\t'+str(currep[0]) + ':' + str(len(orignial[currep_ind][0]))  +'\t'+prestr+ orignial[currep_ind][2]);
                #logging.debug('\t'+str(currep[0]) + ':' + str(len(orignial[currep_ind][0]))  +'\t'+orignial[currep_ind][2]);
                logging.debug('\t'+str(currep[0]) + ' p:' + str(len(currep[1])) +'\t' + prestr+ (currep[1][:aaprinindex]))
        
        #p2, allocr = get2Peaks(rptrue) 

        return [repeatgene, ref_repeat, p2, allocr, len(rptrue), len(rpfalse)+wrongalign]

def getRepeatForKnownGene(gLoc, repeatgene, bamfile, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id):
	repeatgene = repeatgene.lower()
        mgloc = get_gLoc(repeatgene, gLoc);

        gene_start_end = [int(mgloc[1]), int(mgloc[2])]
        repeat_start_end = [int(mgloc[3]), int(mgloc[4])]
        res = (getRepeatForGivenGene(mgloc[0], repeatgene, gene_start_end, repeat_start_end, bamfile, mgloc[5], mgloc[6],isAlign, isupdown, isExtend, unique_file_id, analysis_file_id))
        res.append(mgloc[7])
        return res;

def getRepeat(gLoc, bamfile, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id):
        summary = []

        glkeys = gLoc.keys(); glkeys.sort()
        for glk in glkeys:
                summary.append(getRepeatForKnownGene(gLoc, glk, bamfile, isAlign, isupdown, isExtend, unique_file_id, analysis_file_id));

	return summary;

