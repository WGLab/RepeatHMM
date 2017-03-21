
import re;
import os;
import sys;
import string;
import math;
import copy

import numpy;
import peakutils;

import logging

from scipy.stats import norm

import heapq


def getValues(mdict, minreads):
        valdict = {}
        mdkeys = mdict.keys();
        for mdk in mdkeys:
                if not valdict.has_key(mdict[mdk]):
                        valdict[mdict[mdk]] = 0;
                valdict[mdict[mdk]] += 1;
        valkeys = valdict.keys(); valkeys.sort();
        more = []
        for i in range(1, len(valkeys)+1):
                if valkeys[-i]<minreads: continue;
                if valdict[valkeys[-i]]>=2: more.append(valkeys[-i])
        return more;

def findSameValues(mdict, v, mdkeys):
        start = -1; end = -1;
        for mdk_ind in range(len(mdkeys)-1):
                if mdict[mdkeys[mdk_ind]]==mdict[mdkeys[mdk_ind+1]] and mdict[mdkeys[mdk_ind]]==v:
                        if start==-1: start = mdk_ind
                        end = mdk_ind+1
                if (not mdict[mdkeys[mdk_ind]]==mdict[mdkeys[mdk_ind+1]]) and (not start==-1): break;
        return [start, end]

def reviseDict(mdict, v):
        has_revised = False;
        mdkeys = mdict.keys();  mdkeys.sort()
        while True:
                [start, end] = findSameValues(mdict, v, mdkeys)
                if start==-1 or end==-1: break;
                else:
                        i=start; j = end;
                        while True:
                                if i>=j:
                                        print 'Error i is larger than j', i, j, start, end;
                                        break;
                             
                                has_revised = True;  
                                i += 1;
                                mdict[mdkeys[i]] += 1;
                                if i==j: break;
                                j -= 1;
                                if i==j: break;
                                mdict[mdkeys[j]] += 1;

        return [mdict, has_revised];

def reviseDictAccordingV(mdict, minreads):
	while True:
		more = getValues(mdict, minreads)
		more.sort();
		has_revised = False;
		for m in more:
			mdict, has1_revised = reviseDict(mdict, m);
			if not has_revised: has_revised = has1_revised
		if not has_revised: break;
	return mdict

def getCloseSkipNeighbor(lendict, msign, xi, mindex):
	moreneighers = int(mindex/100.0+0.25);
	reti = None;
	if lendict.has_key(xi):
		reti = xi;
	elif xi>100:
		for fari in range(moreneighers):
			if msign<0 and lendict.has_key(xi-1-fari):
				reti = xi -1-fari ; break;
			elif msign>0 and lendict.has_key(xi+1+fari):
				reti = xi +1+fari; break;
	return reti

def getNeighbors(mindex, lendict):
	neighbor_sum = 0; nei_num = 0;
	moreneighers = int(mindex/100.0+0.25);

	for xi in range(mindex-1, mindex+2):
		reti = getCloseSkipNeighbor(lendict, xi-mindex, xi, mindex)
		if not reti==None:
			neighbor_sum += lendict[reti]; nei_num += 1
		continue;
		
	return neighbor_sum, nei_num

def getNeighbors2(x_index2a, x_index2b, lendict):
	neighbor3 = [0,0, 1, 1]; point3 = [x_index2a, x_index2b]
	for pi in range(len(point3)):
		neighbor3[pi], neighbor3[pi+2] = getNeighbors(point3[pi], lendict)
	
	neighbor3.append(neighbor3[2]-1);
	if neighbor3[4]<1: neighbor3[4]==1
	neighbor3.append(neighbor3[3]-1); 
	if neighbor3[5]<1: neighbor3[5]==1
	
	return neighbor3

def getRatioInfo(neighbor3, x_index2a, x_index2b):
	cur_ratio = round(neighbor3[1]/float(neighbor3[0]),3)
	cur_ratio_threshold = getLargerOverSmallRatio(x_index2a, x_index2b)
	
	return [cur_ratio, cur_ratio_threshold]

def selectFromTwoX(x_index2a, x_index2b, lendict, firstIndex, minreads):
	if x_index2a==x_index2b: x_index2 = x_index2b
	else:
		if x_index2a>x_index2b:
			x_index2a, x_index2b = x_index2b, x_index2a

		neighbor3 = getNeighbors2(x_index2a, x_index2b, lendict)

		msmall1 = False;
		if neighbor3[0]>=minreads*neighbor3[4] or checkSmallSupport1(x_index2a, lendict, minreads): msmall1 = True;
		msmall2 = False;
		if neighbor3[1]>=minreads*neighbor3[5] or checkSmallSupport1(x_index2b, lendict, minreads): msmall2 = True;
		
		#print msmall2, neighbor3, minreads, checkSmallSupport1(x_index2b, lendict, minreads)
		
		closecheck1 = checkClose2Peak(x_index2a, firstIndex, lendict)
		closecheck2 = checkClose2Peak(x_index2b, firstIndex, lendict)
		
		if (x_index2a not in closecheck1) and (x_index2b not in closecheck2): pass
		elif (x_index2a in closecheck1) and (x_index2b in closecheck2): pass
		else:
			if (x_index2a not in closecheck1) and msmall1: msmall1 = False
			if (x_index2b not in closecheck2) and msmall2: msmall2 = False;
		
		#print msmall1, msmall2, x_index2a, x_index2b, closecheck1, closecheck2

		if msmall1 and msmall2:
			ratioInfo = getRatioInfo(neighbor3, x_index2a, x_index2b)
			ratioInfo2= getRatioInfo([lendict[x_index2a], lendict[x_index2b]], x_index2a, x_index2b)

			#print ratioInfo, ratioInfo2

			if ratioInfo[0]<ratioInfo[1]: 
				if ratioInfo2[0]<ratioInfo2[1]:
					x_index2 = x_index2a
				else: x_index2 = x_index2b
			else: x_index2 = x_index2b
		else:
			if msmall1: x_index2 = x_index2a
			elif msmall2: x_index2 = x_index2b
			else:
				if x_index2a<x_index2b:  x_index2 = x_index2b
				else: x_index2 = x_index2a
		
	return x_index2

def selectFromTwo(y_index2a, y_index2b, lendict, y_index1, x, minreads):
	if y_index2a==y_index2b: y_index2 = y_index2b
	else:
		x2 = selectFromTwoX(x[y_index2a], x[y_index2b], lendict, x[y_index1], minreads)
		if x2==x[y_index2a]: y_index2 = y_index2a
		else: y_index2 = y_index2b
	return y_index2
	

def selectOne(iy, curyvalue, lendict, y_index1, indexes, x, minreads):
	y_index2a = indexes[iy.index(curyvalue)];
	iy.reverse();
	y_index2b = indexes[len(iy) - iy.index(curyvalue) - 1]

	if y_index2a==y_index1: return y_index2b
	elif y_index2b==y_index1: return y_index2a
	
	return selectFromTwo(y_index2a, y_index2b, lendict, y_index1, x, minreads)


def getPeaks2(x, y, lendict, mm, mdebug, minreads, smalleratio, issum, len3dict, pnearby_size=3):
	mRepeatTime = 5
	indexes = peakutils.indexes(numpy.array(y), thres=0.0001) - mm
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

	print 'indexes', indexes

	if len(indexes)>1:
		iy = []
		for i in indexes: iy.append(y[i]);
		topNlargest = 1;
		print iy
		while len(iy)>=topNlargest:
			ylargest = heapq.nlargest(topNlargest, iy);
			y_index1 = indexes[iy.index(ylargest[topNlargest-1])];
			#print 'topNlargest1', topNlargest, y_index1, x[y_index1]
			if x[y_index1]>4: break;
			topNlargest += 1
		topNlargest += 1
		has2 = False;
		while len(iy)>=topNlargest:
			ylargest = heapq.nlargest(topNlargest, iy);
			y_index2 = selectOne(copy.deepcopy(iy), ylargest[topNlargest-1], lendict, y_index1, indexes, x, minreads)
			#print 'topNlargest2', topNlargest, y_index2, x[y_index2]
			if x[y_index2]>4: 
				has2 = True;
				break;
			topNlargest += 1
		for i in range(2):
			topNlargest += 1
			while len(iy)>=topNlargest:
				ylargest = heapq.nlargest(topNlargest, iy);
				p_y_index3 = selectOne(copy.deepcopy(iy), ylargest[topNlargest-1], lendict, y_index1, indexes, x, minreads)
				#print 'topNlargest3',topNlargest, p_y_index3, x[p_y_index3]
				if x[p_y_index3]>4: 
					has2 = True;
					y_index2 = selectFromTwo(y_index2, p_y_index3, lendict, y_index1, x, minreads)
					break;

				topNlargest += 1

		#
		#ylargest = heapq.nlargest(2, iy);
		#y_index1 = indexes[iy.index(ylargest[0])];
		#
		#p_y_index2 = selectOne(copy.deepcopy(iy), ylargest[1], lendict, y_index1, indexes, x, minreads)
		#if len(iy)>2:
		#	ylargest = heapq.nlargest(3, iy);
		#	p_y_index3 = selectOne(copy.deepcopy(iy), ylargest[2], lendict, y_index1, indexes, x, minreads)
		#	#print 'p_y_index3', ylargest, p_y_index3, p_y_index2, x[p_y_index3], x[p_y_index2]
		#	y_index2 = selectFromTwo(p_y_index2, p_y_index3, lendict, y_index1, x, minreads)
		#	if len(iy)>3:
		#		ylargest = heapq.nlargest(4, iy);
		#		p_y_index4 = selectOne(copy.deepcopy(iy), ylargest[3], lendict, y_index1, indexes, x, minreads)
		#		y_index2 = selectFromTwo(p_y_index2, p_y_index4, lendict, y_index1, x, minreads)
		#else: y_index2 = p_y_index2
		##print 'y_index2', y_index2, x[y_index2]

		peak2.append(x[y_index1]); 
		for ni in range(x[y_index1]-pnearby_size+1, x[y_index1]+pnearby_size):
			reti = getCloseSkipNeighbor(lendict, ni-x[y_index1], ni, x[y_index1])
			if not reti==None: pnearby.append(reti);
			#pnearby.append(ni);
		if has2:
			peak2.append(x[y_index2])
			for ni in range(x[y_index2]-pnearby_size+1, x[y_index2]+pnearby_size):
				reti = getCloseSkipNeighbor(lendict, ni-x[y_index2], ni, x[y_index2])
				if not reti==None: pnearby.append(reti);
				#pnearby.append(ni);
	else:
		if len(indexes)==1:
			logging.debug("Only find one peak");
			peak2.append(x[indexes[0]])
			for ni in range(x[indexes[0]]-pnearby_size+1, x[indexes[0]]+pnearby_size):
				reti = getCloseSkipNeighbor(lendict, ni-x[indexes[0]], ni, x[indexes[0]])
				if not reti==None: pnearby.append(reti);
				#pnearby.append(ni);
		else:
			logging.debug("Cannot find a peak");


	#so far at most one peak
	if len(peak2)>1 and y[peak2.index(peak2[1])] > [y[peak2.index(peak2[0])]]:
		peak2[0], peak2[1] = peak2[1], peak2[0]
	
	p2len = len(peak2)

	#for two tails
	twotails = [];
	if (x[0] not in pnearby): twotails.append(x[0]);
	if (x[len(x)-1] not in pnearby): twotails.append(x[len(x)-1]);
	
	tailschoosen = []
	for ti in twotails:
		#for two tails
		tinearby = [ti-1, ti+1]; shouldchoose = True;
		if not checkSmallSupport1(ti, lendict, minreads): shouldchoose=False
		else:
			for tnb in tinearby:
				# sepcial checking
				if lendict.has_key(tnb): 
					ratioInfo = getRatioInfo([lendict[tnb], lendict[ti]], tnb, ti)
					if ratioInfo[0]<ratioInfo[1]:
						shouldchoose=False;
		tailschoosen.append(shouldchoose);
		if not shouldchoose: continue;
	
		if not (ti>=mRepeatTime): continue;
			
		for i in range(p2len):
			if selectFromTwoX(ti, peak2[i], lendict, peak2[0], minreads)==ti:
				peak2.insert(i, ti); break;

	if len(peak2)==0:
		for curtt_ind in range(len(twotails)):
			if not tailschoosen[curtt_ind]: continue;
			curtt = twotails[curtt_ind]
			if curtt>=mRepeatTime and checkSmallSupport1(curtt, lendict, minreads): peak2.append(curtt);

 	if mdebug:
		for i in range(len(x)): print ('%d=%d;' % (x[i], lendict[x[i]])),
		print 'mPeak', peak2
		print '#xF=',
		for i in range(len(x)): print ('%3d' % x[i]),
		print ''
		print 'yF =',
		for i in range(len(x)): print ('%3d' % lendict[x[i]]),
		print ''

	if len(peak2)>1: 
		mstr = str(issum)+ '\tIn Peak:'
		if not lendict.has_key(peak2[0]): mstr += (' No %d' % peak2[0]);
		elif not lendict.has_key(peak2[1]): mstr += (' No %d' % peak2[1]);
		else:  mstr += (' ratio=%.3f(%d/%d)' % (lendict[peak2[1]]/float(lendict[peak2[0]]), lendict[peak2[1]], lendict[peak2[0]]))
		#logging.info(mstr)

		if mdebug:
			print issum, ('%.3f=%d/%d: ' % (lendict[peak2[1]]/float(lendict[peak2[0]]), lendict[peak2[1]], lendict[peak2[0]]) ),

		neighbor3 = getNeighbors2(peak2[0], peak2[1], lendict)
		ratioInfo = getRatioInfo(neighbor3, peak2[0], peak2[1])
		ratioInfo2= getRatioInfo([lendict[peak2[0]], lendict[peak2[1]]], peak2[0], peak2[1])

		mstr += ('\n<1=%d, 0=%d> >>> %.3f(%d/%d)' % (peak2[1], peak2[0], neighbor3[1]*peak2[1]/float(neighbor3[0]*peak2[0]), neighbor3[1]*peak2[1], neighbor3[0]*peak2[0]))

		if mdebug:
			print 'ratioInfo', neighbor3, peak2[0], peak2[1], ratioInfo, math.exp(-71/31), math.exp(-peak2[1]/peak2[0])
			print 'ratioInfo2', [lendict[peak2[0]], lendict[peak2[1]]], peak2[0], peak2[1], ratioInfo2, math.exp(-71/31), math.exp(-peak2[1]/peak2[0])
			print mstr
			print neighbor3[1], peak2[1], '=', neighbor3[1]*peak2[1], '/', neighbor3[0], peak2[0], '=', neighbor3[0]*peak2[0], '===', neighbor3[1]*peak2[1]/float(neighbor3[0]*peak2[0]),'ratioInfo=', ratioInfo, ratioInfo2, issum
		logging.info(mstr)

		if ratioInfo[0]<ratioInfo[1] and ratioInfo2[0]<ratioInfo2[1]:
			peak2 = peak2[:1]

	if issum:
		samevalues = []
		ldkeys = lendict.keys(); ldkeys.sort();
		for curp_ind in range(len(peak2)):
			curp = peak2[curp_ind]
			nb3 = [curp]

			moreneighers = int(curp/100.0+0.25)+1;
			for nbni in range(curp-moreneighers, curp+moreneighers+1):
				#print 'nbni', nbni, curp, moreneighers, range(curp-moreneighers, curp+moreneighers+1)
				if nbni in ldkeys: nb3.append(nbni)
			#pin = ldkeys.index(curp)
			#if pin-1>0: nb3.append(ldkeys[pin-1])
			#if pin+1<len(ldkeys): nb3.append(ldkeys[pin+1])
			#print 'nb3', nb3
			curmx = nb3[-1]
			for curi in nb3:
				ratioInfo = getRatioInfo([lendict[curi], lendict[curmx]], curi, curmx)
				if mdebug: print 'in issum: old', curmx, '>>', curi, curmx, [lendict[curi], lendict[curmx]], ratioInfo, 
				
				checkcuri = checkSmallSupport1(curi, lendict, minreads)
				checkcurmx = checkSmallSupport1(curmx, lendict, minreads)
				if mdebug: print checkcuri, checkcurmx,
				if checkcuri and checkcurmx:
					if ratioInfo[0]<ratioInfo[1]:
						curmx = curi
				elif checkcuri and (not checkcurmx):
					curmx = curi

				if mdebug: print 'new ', curmx

				# special checking
				#if (lendict[curi]*curi>lendict[curmx]*curmx): curmx = curi
				#if (lendict[curi]*curi==lendict[curmx]*curmx) and curi>curmx: curmx = curi
				#if lendict[curi]>lendict[curmx]: curmx = curi
				#elif lendict[curi]==lendict[curmx] and curi>curmx: curmx = curi
			peak2[curp_ind] = curmx

			samevalues.append([curmx, nb3])

	if mdebug: print '1>>', peak2, minreads
	peak2 = checkSmallSupport(peak2, lendict, minreads)
	if mdebug: print '2>>', peak2

	peak2.sort()
	if len(peak2)>1:
		if mdebug:
			print 'compare2', len3dict[peak2[0]], len3dict[peak2[1]], smalleratio, len3dict[peak2[1]]*smalleratio, ' or ', lendict[peak2[0]], lendict[peak2[1]], smalleratio, lendict[peak2[1]]*smalleratio, 'peak2', peak2
		if len3dict[peak2[0]]<len3dict[peak2[1]]*smalleratio or lendict[peak2[0]]<lendict[peak2[1]]*smalleratio:
			peak2 = peak2[1:]

	if len(peak2)>1 and peak2[0]==peak2[1]: peak2=peak2[1:]
	if len(peak2)>2: peak2 = peak2[:2];

	peak2 = checkClosePeak(peak2, lendict)

	if issum:
		return peak2, samevalues
	else: return peak2

def checkClose2Peak(p1, p2, lendict):
	istest = True; istest = False; 
	peak2 = [p1, p2]
	if 1<abs(peak2[0]-peak2[1])<=3 and (-1 not in peak2):
		closep = [peak2[0], peak2[1]]; closep.sort();
		avg2 = 0; minsup = -1;
		for betp in range(closep[0]+1, closep[1]):
			if lendict.has_key(betp): cursup = lendict[betp]
			else: cursup = 0;
			avg2 += cursup;
			if (minsup==-1 or minsup>cursup): minsup = cursup
		avg2 /= 2.0; #
		
		if abs(peak2[0]-peak2[1])==2:
			mycriterion = [avg2,'avg'];
			neigh_perc = 0.45;
			thr0 = lendict[closep[0]]*neigh_perc; thr1 = lendict[closep[1]]*neigh_perc;
		else:
			if minsup==-1: print 'Warning!!!! minimum support between is -1', peak2
			mycriterion = [minsup,'min']
			neigh_perc = 0.8;
			thr0 = lendict[closep[0]]*neigh_perc; thr1 = lendict[closep[1]]*neigh_perc;
		if istest: print 'checkClose2Peak:', mycriterion, thr1, thr0, closep[1], closep[0]
		if mycriterion[0]>=thr1 and mycriterion[0]>=thr0: 
			ratioInfo = getRatioInfo([lendict[closep[0]], lendict[closep[1]]], closep[0], closep[1])
			if istest: print 'ratioInfo=',ratioInfo, [lendict[closep[0]], lendict[closep[1]]], closep[0], closep[1]
			if ratioInfo[0]<ratioInfo[1]:
				return [closep[0]]
			else: return [closep[1]]
		elif mycriterion[0]>=thr1: return [closep[0]]
		elif mycriterion[0]>=thr0: return [closep[1]]
		else: return peak2
	else: return peak2
def checkClosePeak(peak2, lendict):
	newpeak = [];
	if len(peak2)>1:
		newpeak = checkClose2Peak(peak2[0], peak2[1], lendict)
		if len(peak2)>2:
			newpeak.extend(peak2[2:])
	else: newpeak = peak2

	return newpeak

def getLargerOverSmallRatio(p1, p2):
	mp = [int(p1),int(p2)]; #

	#if math.fabs(mp[0]-mp[1])>=3:
	if math.fabs(mp[0]-mp[1])>=5:
		mratio = 3**(-mp[1]/float(mp[0]))
	else:
		if abs(mp[0]-mp[1])==4: mratio=0.65;
		elif abs(mp[0]-mp[1])==3: mratio=0.75;
		elif abs(mp[0]-mp[1])==2: mratio=0.85; #mratio=0.7
		else: mratio=0.95
		#mpdif = mp[0]-mp[1]
		#if mpdif==-2: mratio=0.55
		#elif mpdif==2: mratio=0.7
		#elif mpdif==-1: mratio=0.85
		#else: mratio=0.95

	if mratio<0.001: mratio=0.001;

	return mratio

def checkSmallSupport1(msup, lendict, minreads):
	#print 'checkSmallSupport1', msup, minreads, lendict[msup]
	if lendict[msup]<minreads:
		return False;
	else: return True;

def checkSmallSupport(peak2, lendict, minreads):
	newpeak2 = []
	for p in peak2:
		if checkSmallSupport1(p, lendict, minreads):
			newpeak2.append(p)
	return newpeak2

def get2Peaks(lengd, MinSup, mdebug = False, commonoptions=None):
	lendict = {};
	for l in lengd:
		l = int(l+0.5)
		if not lendict.has_key(l): lendict[l] = 0;
		lendict[l] += 1;
	
	return get2PeaksfromDict(lendict, MinSup, mdebug, commonoptions)

def getNewRepeatForIllumina(lendict, MinSup, commonoptions, oldp2):
   if not ( len(oldp2)==0 or (len(oldp2)==1 and (not oldp2[0]==0)) or (len(oldp2)>1 and oldp2[0]==oldp2[1] and (not oldp2[0]==0)) ): return oldp2
   
   newrepeatsdictkeys = lendict.keys(); newrepeatsdictkeys.sort()
   max1 = [0,0]; max2 = [0,0]
   allreads = 0;
   for nrkey in newrepeatsdictkeys:
      allreads = allreads + lendict[nrkey]
      if lendict[nrkey]>=MinSup:
         if lendict[nrkey]>max1[1]:
            max2 = [max1[0], max1[1]]
            max1 = [nrkey, lendict[nrkey]]
         else:
            if lendict[nrkey]>max2[1]:
               max2 = [nrkey, lendict[nrkey]]
   if max1[1]<MinSup:
      print 'Warning!!! less than '+str(MinSup), max1, max2
   print max1, max2, max2[1]/float(max1[1]), (max2[1]+max1[1])/float(allreads)
   if max2[1]==0:
      newrepeats = [max1[0], max1[0]]
   else:
      if max2[1]/float(max1[1])>0.6 or (max2[1]+max1[1])/float(allreads)>0.95:
         newrepeats = [max1[0], max2[0]]
      else: newrepeats = [max1[0], max1[0]]

   newrepeats.sort()

   if len(oldp2)>0 and (oldp2[0] not in newrepeats):
      print 'Warning!!! not in ', oldp2, newrepeats
 
   return newrepeats


def get2PeaksfromDict(lendict, MinSup, mdebug = False, commonoptions=None):
	isillumina = False;
	if (not commonoptions==None) and commonoptions['SeqTech']=="Illumina":
		isillumina = True;

	minreads = int(MinSup);
	if minreads<2: minreads = 2

	smalleratio = 0.5; #0.55 #0.4; #0.55 # 0.5 #0.6;

	
	ldkeys = lendict.keys(); ldkeys.sort();
	
	allocr = 'allocr:';
	for ldk in ldkeys:
		allocr += ('%d:%d, ' % (ldk, lendict[ldk]))
	if mdebug: print allocr
	logging.info(allocr)
	print allocr, minreads, MinSup
	
	for ldk in ldkeys:
		if ldk<4: del lendict[ldk]
	ldkeys = lendict.keys(); ldkeys.sort();
	
	#for special case;
	if len(ldkeys)<1: return [[0], allocr]
	elif len(ldkeys)<2: peak2 = [ldkeys[0]]; #[[ldkeys[0]], allocr]; 
	elif len(ldkeys)==2 or len(ldkeys)==3:
		maxk = ldkeys[0]
		if maxk==0: maxk = ldkeys[1]
		for ldk in ldkeys:
			if lendict[ldk]>=lendict[maxk] and ldk>0: maxk = ldk
		peak2 = [maxk]
	if len(ldkeys)<4:
		peak2 = checkSmallSupport(peak2, lendict, minreads)
		if len(peak2)==0: peak2 = [0]

		if isillumina:
			peak2 = getNewRepeatForIllumina(lendict, MinSup, commonoptions, peak2)

		return [peak2, allocr[:-1]];

	#get robustness list with a window size=3;
	len3dict = {}
	for i in range(len(ldkeys)):
		curlk = ldkeys[i]
		len3dict[curlk] = getNeighbors(curlk, lendict)[0]
		continue;
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
	
	len3dict = reviseDictAccordingV(len3dict, minreads)
	
	x = []; yo = [];
	for ldk in ldkeys:
		x.append(ldk); 
		yo.append(len3dict[ldk]);
	peak2, samevalues = getPeaks2(x, yo, lendict, 0, mdebug, minreads, smalleratio, True, len3dict)
	if mdebug: print peak2

	peaksteststr = 'Peak2='+str(len(peak2))+' info:';
	for pst in peak2: peaksteststr += ' '+str(pst)
	logging.info(peaksteststr)
	
	if len(peak2)==1:
		curPeak1 = peak2[0]
		lenP1dict = {}
		for px in range(curPeak1-6, curPeak1+7):
			if lendict.has_key(px):
				lenP1dict[px] = lendict[px];
			else: pass 
		ld1keys = lenP1dict.keys(); ld1keys.sort();
		#
		forbelowtest = False; #True;
		if forbelowtest: #########################comment
			print 'test', int(minreads/2.0+0.5), 
			allocr = '';
			for ldk in ld1keys:
				allocr += ('%d:%d, ' % (ldk, lenP1dict[ldk]))
			print allocr; 
		lenP1dict = reviseDictAccordingV(lenP1dict, int(minreads/2.0+0.5))
		#
		if forbelowtest: #########################comment
			print 'test', int(minreads/2.0+0.5),
			allocr = '';
			for ldk in ld1keys:
				allocr += ('%d:%d, ' % (ldk, lenP1dict[ldk]))
			print allocr
		
		ld1keys = lenP1dict.keys(); ld1keys.sort();
		x = []; yo = []
		for ldk in ld1keys:
			x.append(ldk);
			yo.append(lenP1dict[ldk]);
		peak2 = getPeaks2(x, yo, lendict, 0, mdebug, int(minreads/2.0+0.5), smalleratio, False, len3dict, 0)

		peaksteststr = 'Peak2 for P1='+str(len(peak2))+' info:';
		for pst in peak2: peaksteststr += ' '+str(pst)
		logging.info(peaksteststr)
	
		if curPeak1 not in peak2:
			allocr = '';
			for ldk in ld1keys:
				allocr += ('%d:%d, ' % (ldk, lenP1dict[ldk]))
			
			p1pstr = ' '
			for p1p in peak2: p1pstr += str(p1p)+','
			logging.info("Warning!!!! previously detected peak "+str(curPeak1)+" not in this detection:"+p1pstr[1:-1]+" for "+allocr)
			
			newpeak2 = []; 
			for curdetp in peak2:
				issame = False;
				for i in range(len(samevalues)):
					if (curPeak1 in samevalues[i][1]) and (curdetp in samevalues[i][1]) and (lendict[curPeak1]==lendict[curdetp]):
						issame = True; newpeak2.append(samevalues[i][0])
				if not issame: newpeak2.append(curdetp)
			peak2 = newpeak2
			p1pstr = ' '
			for p1p in peak2: p1pstr += str(p1p)+','
			logging.info("New peaks "+p1pstr[1:-1])
	
	if len(peak2)>1 and peak2[0]==peak2[1]: peak2=peak2[1:]
	if len(peak2)>2: peak2 = peak2[:2]
	if len(peak2)<2:
		if lendict[ldkeys[-1]]>minreads and (len(peak2)==0 or (len(peak2)==1 and math.fabs(peak2[0]-ldkeys[-1])>1 and ((not lendict.has_key(ldkeys[-1]-1)) or (lendict.has_key(ldkeys[-1]-1) and lendict[ldkeys[-1]]>lendict[ldkeys[-1]-1]) ) )):
			peak2.append(ldkeys[-1])

	peak2.sort();
	if len(peak2)>1:
		neighbor3 = getNeighbors2(peak2[0], peak2[1], lendict)
		ratioInfo = getRatioInfo(neighbor3, peak2[0], peak2[1])
		ratioInfo2= getRatioInfo([lendict[peak2[0]], lendict[peak2[1]]], peak2[0], peak2[1])

		mstr = ('P1 <1=%d, 0=%d> >>> %.3f(%d/%d)' % (peak2[1], peak2[0], neighbor3[1]*peak2[1]/float(neighbor3[0]*peak2[0]), neighbor3[1]*peak2[1], neighbor3[0]*peak2[0]))
		logging.info(mstr)

		if mdebug: 
			print mstr
			print neighbor3, ratioInfo

		if ratioInfo[0]<ratioInfo[1] and ratioInfo2[0]<ratioInfo2[1]:
			peak2 = peak2[:1]
	
	if mdebug: print 'mPeak', peak2

	peak2 = checkSmallSupport(peak2, lendict, minreads)
	peak2 = checkClosePeak(peak2, lendict)

	if isillumina:
		peak2 = getNewRepeatForIllumina(lendict, MinSup, commonoptions, peak2)

	return [peak2, allocr[:-1]];

if __name__=='__main__':
	mydict1 = {12: 6, 13: 9, 14: 41, 15: 79, 16: 55, 17: 16, 18: 4, 19: 2, 21: 1, 23: 4, 26: 2, 28: 1, 29: 1, 30: 1, 31: 1, 32: 2, 33: 2, 34: 3, 35: 1, 36: 2, 40: 3, 41: 1, 42: 1, 43: 2, 45: 2, 47: 3, 49: 2, 50: 2, 51: 2, 52: 1, 54: 1, 55: 1, 56: 3, 57: 1, 58: 1, 59: 2, 61: 1, 63: 2, 68: 2, 69: 1, 71: 1, 73: 1, 75: 1, 77: 1, 78: 1, 79: 2, 80: 7, 81: 28, 82: 66, 83: 49, 84: 50, 85: 41, 86: 28, 87: 6, 88: 1}

	mydict2 = {22:4, 23:5, 24:16, 25:55, 26:67, 27:23, 28:6, 29:12, 30:4, 31:4, 32:2, 33:2, 43:1, 48:2, 50:1, 53:1, 55:1, 56:1, 59:1, 60:1, 61:1, 62:3, 63:3, 64:1, 65:1, 66:3, 67:3, 68:6, 69:12, 70:20, 71:35, 72:35, 73:35, 74:15, 75:6, 76:6, 77:2, 78:4}

	mydict3 = {39:1, 40:3, 41:4, 42:4, 43:4, 44:1, 45:1, 46:1, 47:1, 56:1, 61:1, 71:1, 72:2, 73:2, 74:4, 75:4, 76:4, 78:1}

	mydict4 = {46:1, 47:1, 49:1, 50:1, 51:9, 52:11, 53:18, 54:19, 55:10, 56:7, 57:3, 58:3, 59:1, 60:2, 61:1, 62:1, 63:3, 64:4, 65:5, 66:17, 67:10, 68:19, 69:5, 70:3, 71:2, 72:2, 73:1}

	mydict5 = {58:2, 59:2, 61:5, 62:3, 63:6, 64:1, 65:1, 68:1, 75:2, 76:4, 77:3, 79:5, 80:1, 82:2, 83:2}

	#print get2PeaksfromDict(mydict4, 2);

	mydict6 = {7:1, 9:1, 10:1, 11:1, 12:1, 13:9, 14:6, 15:8, 16:7, 17:13, 18:10, 19:7, 20:1}
	#print get2PeaksfromDict(mydict6, 10), '\n'

	mydict7 = {0:1, 4:1, 11:1, 18:1, 19:3, 20:1, 24:2, 27:3, 28:9, 29:5, 30:4, 31:6, 32:6, 33:6, 34:3, 35:5, 36:2, 37:4,}
	#print get2PeaksfromDict(mydict7, 10), '\n';
	
	mydict8 = {16:2, 17:1, 19:4, 20:2, 21:6, 22:10, 23:10, 24:10, 25:19, 26:6, 27:6, 28:2, 31:1}
	#print get2PeaksfromDict(mydict8, 10), '\n';

	mydict9 = {19:1, 30:1, 36:1, 37:1, 38:2, 40:2, 41:1, 42:3, 43:1, 44:8, 45:3, 46:8, 47:9, 48:8, 49:7, 50:2, 51:5, 52:2, 53:5, 54:1, 55:2, 57:2, 60:1, 62:2, 66:1}
	#print get2PeaksfromDict(mydict9, 10), '\n';

	mydict10 = {5:1, 6:3, 7:11, 8:24, 9:29, 10:17, 11:28, 12:1}
	mydict10 = {18:1, 21:1, 22:1, 23:2, 24:1, 25:6, 26:5, 27:2, 28:6, 29:2, 30:4, 31:5, 32:1, 33:5, 34:6, 35:6, 36:3, 37:6, 38:3, 39:1, 73:1}
	mydict10 = {0:1, 8:1, 9:3, 10:1, 13:5, 14:5, 15:5, 16:12, 17:25, 18:20, 19:14, 20:5}
	mydict10 = {7:1, 8:6, 9:22, 10:17, 11:27, 12:26, 13:1, 14:1}
	mydict10 = {10:8, 11:4, 12:9, 13:11, 14:18, 15:14, 16:17, 17:6, 18:1}
	mydict10 = {18:1, 21:1, 22:1, 23:2, 24:1, 25:6, 26:5, 27:2, 28:6, 29:2, 30:4, 31:5, 32:1, 33:5, 34:6, 35:6, 36:3, 37:6, 38:3, 39:1, 73:1}
	mydict10 = {14:1, 23:1, 24:3, 25:3, 27:3, 28:9, 29:6, 30:15, 31:7, 32:20, 33:11, 34:4, 35:7, 36:2, 37:3, 38:1, 39:2, 40:1}
	mydict10 = {7:1, 9:1, 10:1, 11:1, 12:1, 13:9, 14:6, 15:8, 16:7, 17:13, 18:10, 19:7, 20:1}
	mydict10 = {6:1, 7:2, 8:4, 9:8, 10:11, 11:13, 12:10, 13:5, 14:11, 15:16, 16:11, 17:5, 18:2, 19:2}
	mydict10 = {0:1, 11:3, 13:1, 14:1, 15:13, 16:11, 17:10, 18:10, 19:6, 20:6, 21:12, 22:7, 23:2, 24:1, 25:1}
	mydict10 = {8:2, 9:7, 10:8, 11:7, 12:12, 13:8, 14:9, 15:6, 17:2}
	mydict10 = {0:2, 7:1, 9:2, 10:4, 11:10, 12:9, 13:8, 14:10, 15:9, 16:13, 17:8, 18:12, 19:3, 20:4, 21:4, 22:3, 23:3, 24:5, 26:1}
	mydict10 = {10:4, 11:10, 12:11, 13:6, 14:12, 15:10, 16:8, 17:4, 18:2}
	mydict10 = {10:4, 11:9, 12:11, 13:8, 14:13, 15:17, 16:15, 17:1}
	mydict10 = {0:2, 7:1, 9:2, 10:4, 11:10, 12:9, 13:8, 14:10, 15:9, 16:13, 17:8, 18:12, 19:3, 20:4, 21:4, 22:3, 23:3, 24:5, 26:1}
	mydict10 = {18:1, 21:1, 22:1, 23:2, 24:1, 25:6, 26:5, 27:2, 28:6, 29:2, 30:4, 31:5, 32:1, 33:5, 34:6, 35:6, 36:3, 37:6, 38:3, 39:1, 73:1}
	mydict10 = {6:1, 13:2, 16:1, 17:4, 18:1, 19:8, 20:3, 21:1, 22:1, 23:1, 25:1, 26:3, 27:4, 28:2, 29:4, 30:4, 31:1}
	mydict10 = {12:1, 15:1, 25:1, 27:1, 28:15, 29:88, 30:622, 31:622, 32:112, 33:10, 35:3, 36:2, 37:1, 38:5, 39:2, 41:4, 42:7, 43:3, 44:6, 45:5, 46:7, 47:9, 48:14, 49:5, 50:6, 51:9, 52:7, 53:14, 54:9, 55:8, 56:8, 57:12, 58:7, 59:9, 60:13, 61:6, 62:2, 63:8, 64:10, 65:1, 66:5, 67:15, 68:12, 69:31, 70:50, 71:54, 72:31, 73:19, 74:10, 75:5, 76:3, 77:2, 84:1}
	mydict10 = {3:1, 4:3, 8:1, 10:3, 11:5, 12:7, 13:8, 14:11, 15:8, 16:8, 17:5, 18:4, 19:2, 20:1, 21:7, 22:8, 23:6, 24:5, 25:6, 26:5, 27:3}
	mydict10 = {0:37, 1:3, 2:9, 3:17, 4:14, 5:16, 6:32, 7:34, 8:37, 9:39, 10:41, 11:55, 12:62, 13:94, 14:118, 15:180, 16:513, 17:837, 18:180, 19:56, 20:4, 24:1, 28:1, 31:1, 52:1, 88:1, 104:1, 159:1, 202:1, 204:1, 229:1, 262:1, 284:1, 311:1, 316:1, 326:1, 337:1, 340:1, 344:1, 388:1, 390:1, 417:1, 440:1, 452:1, 468:1, 470:1, 481:1, 485:1, 492:1, 494:1, 495:1, 496:1, 501:1, 521:2, 530:1, 533:1, 534:1, 544:1, 546:1, 549:1, 550:1, 554:1, 558:2, 561:1, 562:1, 563:2, 564:1, 572:1, 574:2, 578:1, 580:1, 588:1, 596:1, 598:1, 599:1, 600:2, 602:1, 603:1, 606:1, 607:2, 608:2, 611:2, 613:2, 614:2, 615:2, 616:2, 619:4, 623:1, 624:1, 625:1, 627:2, 628:1, 629:2, 635:1, 637:1, 638:1, 639:2, 641:2, 642:1, 643:2, 644:1, 645:2, 646:2, 647:1, 648:2, 649:1, 652:2, 653:1, 655:2, 656:2, 657:1, 658:3, 659:4, 660:3, 661:2, 662:1, 663:2, 666:1, 667:2, 668:2, 669:1, 670:3, 672:1, 673:3, 674:1, 675:4, 676:1, 677:3, 678:2, 680:3, 681:2, 682:1, 683:4, 685:2, 686:3, 687:2, 688:1, 690:2, 691:2, 692:1, 693:1, 694:1, 695:3, 696:1, 697:4, 698:3, 700:1, 702:2, 703:4, 704:2, 705:2, 706:2, 707:1, 708:2, 709:4, 710:2, 712:1, 713:5, 715:2, 716:1, 717:2, 718:6, 719:2, 720:2, 722:5, 723:6, 724:1, 725:1, 726:1, 727:3, 728:3, 730:1, 731:5, 732:2, 733:6, 734:2, 736:4, 737:1, 738:2, 739:2, 740:1, 741:4, 742:2, 743:1, 744:5, 745:4, 746:5, 747:4, 748:4, 749:1, 750:2, 751:3, 752:4, 753:3, 754:2, 755:1, 756:4, 757:3, 758:1, 759:2, 760:3, 761:3, 762:3, 764:3, 765:6, 766:3, 767:6, 768:8, 769:4, 770:5, 771:2, 772:9, 773:1, 774:2, 775:6, 776:4, 777:3, 778:4, 779:2, 780:4, 781:5, 782:3, 783:6, 784:4, 785:7, 786:4, 787:3, 788:4, 789:4, 790:7, 791:4, 792:5, 793:4, 794:7, 795:7, 796:5, 797:5, 798:6, 799:3, 800:3, 801:5, 802:5, 803:4, 804:4, 805:6, 806:5, 807:5, 808:3, 809:1, 810:3, 811:1, 812:1, 813:2, 814:1, 815:1, 816:1, 817:1, 818:3, 819:1, 821:1, 822:1, 824:1, 825:1, 836:1, 865:1, 876:1, 914:1, 926:1}
	#print get2PeaksfromDict(mydict10, 10, True), '\n';
	#print get2PeaksfromDict(mydict10, 2, True), '\n';

	mydict10 = {0:1, 2:1, 4:2, 5:6, 6:4, 7:20, 8:29, 9:10, 10:19, 11:11, 12:2}
	mydict10 = {5:2, 6:4, 7:18, 8:1, 9:6, 10:3, 11:20, 12:88, 13:178, 14:38, 15:3, 113:1, 114:1, 117:2, 118:4, 119:7, 120:9, 121:7, 122:3, 123:3, 124:1}
	mydict10 = {5:1, 6:5, 7:6, 8:11, 9:4, 10:6, 11:32, 12:81, 13:33, 14:3, 15:1, 125:2, 126:1, 127:2, 130:5, 132:4, 133:1, 134:1}
	mydict10 = {5:2, 6:4, 7:18, 8:1, 9:6, 10:3, 11:20, 12:88, 13:178, 14:38, 15:3, 113:1, 114:1, 117:2, 118:4, 119:7, 120:9, 121:7, 122:3, 123:3, 124:1}
	mydict10 = {22:1, 23:1, 24:4, 25:3, 26:5, 27:2, 28:13, 29:50, 30:50, 31:22, 32:4, 33:2, 102:1, 103:2, 104:3, 105:3, 106:5, 107:4, 108:9, 109:8, 110:5, 112:2}
	mydict10 = {13:1, 14:1, 15:6, 16:5, 17:7, 18:1, 19:18, 20:48, 21:67, 22:14, 23:5, 126:3, 128:3, 130:1, 132:5, 133:2, 134:2, 135:4, 136:2, 137:2, 138:2}
	mydict10 = {0:5, 3:11, 4:47, 5:203, 6:326, 7:47, 8:21, 9:7, 10:4, 11:2, 12:6, 13:2, 15:1, 222:1, 224:1, 226:1, 227:5, 228:1, 229:3, 230:2, 231:2, 232:1}
	mydict10 = {0:11, 3:15, 4:69, 5:292, 6:455, 7:71, 8:28, 9:11, 10:5, 11:11, 12:2, 13:3, 216:2, 218:7, 219:3, 220:3, 221:6, 222:2, 223:2, 227:1}
	mydict10 = {0:32677, 1:187, 2:2212, 3:837, 4:438, 5:249, 6:173, 7:106, 8:98, 9:78, 10:78, 11:66, 12:61, 13:54, 14:65, 15:56, 16:68, 17:51, 18:74, 19:75, 20:93, 21:117, 22:206, 23:343, 24:360, 25:293, 26:179, 27:126, 28:114, 29:91, 30:88, 31:86, 32:86, 33:92, 34:71, 35:69, 36:90, 37:91, 38:103, 39:95, 40:98, 41:76, 42:73, 43:45, 44:42, 45:25, 46:8, 47:11, 48:7, 49:3, 50:1, 51:1, 54:1}
	mydict10 = {10:4, 11:13, 12:70, 13:310, 14:1132, 15:2551, 16:339, 17:105, 18:44, 19:46, 20:39, 21:26, 22:9, 218:1, 219:3, 220:3, 221:2, 222:5, 223:1, 224:1, 226:1}
	mydict10 = {0:10, 3:26, 4:87, 5:367, 6:1305, 7:2289, 8:336, 9:144, 10:57, 11:30, 12:20, 13:22, 14:13, 15:2, 16:1, 17:1, 108:1, 111:3, 112:2, 113:1, 114:1, 115:1, 116:1}
	mydict10 = {10:1, 11:18, 12:84, 13:318, 14:1121, 15:2534, 16:350, 17:115, 18:44, 19:38, 20:57, 21:23, 22:5, 23:4, 24:1, 240:1, 241:3, 242:2, 243:2, 245:1}
	mydict10 = {8:1, 9:1, 10:5, 11:7, 12:21, 13:78, 14:236, 15:651, 16:1302, 17:1889, 18:273, 19:81, 20:44, 21:44, 22:34, 23:13, 243:2, 244:6, 245:4, 246:4, 247:1, 248:1, 249:2}
	mydict10 = {0:212, 3:98, 4:203, 5:332, 6:760, 7:1393, 8:2017, 9:948, 10:1407, 11:1379, 12:3282, 13:6315, 14:6279, 15:1477, 16:244, 17:77, 18:34, 19:32, 20:24, 21:12, 22:5, 23:4, 24:5, 25:2, 26:1, 27:1, 28:2, 29:2, 31:6, 32:5, 33:3, 34:1, 35:5, 37:5, 38:3, 39:1, 40:1, 48:1, 50:1, 51:1, 55:2, 59:1, 60:1, 61:1, 62:1, 67:3, 68:3, 69:1, 70:3, 71:2, 72:2, 73:4, 74:4, 75:8, 76:6, 77:3, 78:6, 79:9, 80:13, 81:12, 82:18, 83:40, 84:48, 85:61, 86:51, 87:73, 88:72, 89:98, 90:119, 91:98, 92:97, 93:69, 94:65, 95:44, 96:37, 97:20, 98:15, 99:6, 100:5, 101:5, 102:3, 103:2, 104:1, 105:1, 106:1, 107:1, 108:1, 111:1, 113:1, 120:1, 132:1, 152:1, 185:1}
	mydict10 = {7:60, 9:49}; commonoptions = {}; commonoptions['SeqTech']="Illumina"
	mydict10 = {5:23, 13:52}
	mydict10 = {8:1, 13:47, 15:48, 16:4}
	mydict10 = {16:1, 31:1, 34:12, 35:9}
	mydict10 = {0:102, 1:2, 2:2, 3:4, 4:11, 5:6, 6:3, 7:7, 8:3, 9:6, 10:11, 11:8, 12:7, 13:4, 14:6, 15:5, 16:3, 17:3, 20:1, 22:1, 23:1}
	print get2PeaksfromDict(mydict10, 5, True, commonoptions), '\n';

	sys.exit()
	
	dlist = [mydict1, mydict2, mydict3, mydict4, mydict5]
	for dl in dlist:
		dlkeys = dl.keys(); dlkeys.sort();
		for dlk in dlkeys:
			print ('%d:%d' % (dlk, dl[dlk])),
		print ''
		for cov in range(10, 110, 10):
			print '\tcov=', cov, get2PeaksfromDict(dl, cov/10)[0]


