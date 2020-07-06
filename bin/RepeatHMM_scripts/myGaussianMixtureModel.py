
import os;
import sys;
import string;
import math;

import logging

import numpy as np
from sklearn.mixture import GaussianMixture as GM #GMM

#from .myheader import *
from . import myheader

def getLargerOverSmallRatio(p1, p2):
	mp = [int(p1),int(p2)]; #

	if math.fabs(mp[0]-mp[1])>=5:
		mratio = 3**(-mp[1]/float(mp[0]))
	else:
		if abs(mp[0]-mp[1])==4: mratio=0.50;
		elif abs(mp[0]-mp[1])==3: mratio=0.65;
		elif abs(mp[0]-mp[1])==2: mratio=0.80;
		else: mratio=0.95

	if mratio<0.001: mratio=0.001;

	return mratio

def getRatioInfo(x_index2a, x_index2b, a_sum, b_sum):
	cur_ratio = round(b_sum/float(a_sum),3)
	cur_ratio_threshold = getLargerOverSmallRatio(x_index2a, x_index2b)

	return [cur_ratio, cur_ratio_threshold]

def getWindowForCounts(cur_point):
	cur_window = int(cur_point/200.0+0.75) + 1
	return cur_window

def getNeighbors_reads_fixed(lendict, cur_point, cur_window):
	cur_sum_read = 0; total_k = cur_window*2+1; has_k = 0;
	for curk in range(cur_point-cur_window, cur_point+cur_window+1):
		if lendict.has_key(curk):
			cur_sum_read += lendict[curk]
			has_k += 1;
	return [cur_sum_read, has_k, total_k]

def getNeighbors_reads(lendict, cur_point, minreads):
	cur_window = getWindowForCounts(cur_point) #int(cur_point/200.0+0.75)
	cur_sum_read = 0; total_k = cur_window*2+1; has_k = 0;
	for curk in range(cur_point-cur_window, cur_point+cur_window+1):
		if lendict.has_key(curk):
			cur_sum_read += lendict[curk]
			has_k += 1;
	#print cur_point, [cur_sum_read, has_k, total_k]
	
	#if cur_window==0:
	#	for i in [-1, 1]:
	#		if lendict.has_key(cur_point+i): 
	#			print cur_point+i, lendict[cur_point+i], cur_sum_read
	#			if lendict[cur_point+i]>lendict[cur_point]: 
	#				if lendict[cur_point+i]<cur_sum_read:
	#					print 'Error!!!! ', lendict[cur_point+i], cur_sum_read, cur_point, i
	#				cur_sum_read = lendict[cur_point+i]
	
	return [cur_sum_read, has_k, total_k]

def checkSmallSupport1(msup, lendict, minreads):
	cursum = getNeighbors_reads(lendict, msup, minreads)
	if cursum[0]/2<minreads:
		return False;
	else: return True;

def checkSmallSupport(peak2, lendict, minreads):
	newpeak2 = []
	for p in peak2:
		if checkSmallSupport1(p, lendict, minreads):
			newpeak2.append(p)
	return newpeak2

def getNeighbors2(x_index2a, x_index2b, lendict, minreads):
	neighbor3 = [0,0, 1,1,1,1]; point3 = [x_index2a, x_index2b]
	for pi in range(len(point3)):
		neighbor3[pi], neighbor3[pi+2], neighbor3[pi+4] = getNeighbors_reads(lendict, point3[pi], minreads)

	return neighbor3

def getClose(lendict, cur_point, minreads):
   cur_window = getWindowForCounts(cur_point) #int(cur_point/200.0+0.75)
   cur_max_red = 0;
   for curk in range(cur_point-cur_window, cur_point+cur_window+1):
      if lendict.has_key(curk):
         if lendict[curk]>cur_max_red: cur_max_red = lendict[curk]
   return cur_max_red

def selectFromTwoX(x_index2a, x_index2b, lendict, minreads):
	if x_index2a>x_index2b:
		x_index2a, x_index2b = x_index2b, x_index2a

	neighbor3 = getNeighbors2(x_index2a, x_index2b, lendict, minreads)

	msmall1 = False;
	if checkSmallSupport1(x_index2a, lendict, minreads): msmall1 = True;
	msmall2 = False;
	if checkSmallSupport1(x_index2b, lendict, minreads): msmall2 = True;

	if msmall1 and msmall2: #getRatioInfo(x_index2a, x_index2b, a_sum, b_sum):
		ratioInfo = getRatioInfo(x_index2a, x_index2b, neighbor3[0], neighbor3[1])
		risingle = getRatioInfo(x_index2a, x_index2b, getClose(lendict, x_index2a, minreads), getClose(lendict, x_index2b, minreads))
		#print ratioInfo, neighbor3
		#print risingle, getClose(lendict, x_index2a, minreads), getClose(lendict, x_index2b, minreads)
		if ratioInfo[0]<ratioInfo[1] and risingle[0]<risingle[1]:
		#if ratioInfo[0]<ratioInfo[1]:
			x_index2 = x_index2a
		else: x_index2 = x_index2b
	else:
		if msmall1: x_index2 = x_index2a
		elif msmall2: x_index2 = x_index2b
		else:
			if x_index2a<x_index2b:  x_index2 = x_index2b
			else: x_index2 = x_index2a

	return x_index2


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
		if commonoptions['outlog'] <= myheader.M_WARNING: print ('Warning!!! less than '+str(MinSup), max1, max2)
	if commonOptions['outlog'] <= myheader.M_INFO: print (max1, max2, max2[1]/float(max1[1]), (max2[1]+max1[1])/float(allreads))
	if max2[1]==0:
		newrepeats = [max1[0], max1[0]]
	else:
		if max2[1]/float(max1[1])>0.6 or (max2[1]+max1[1])/float(allreads)>0.8:
			newrepeats = [max1[0], max2[0]]
		else: newrepeats = [max1[0], max1[0]]

	newrepeats.sort()

	if len(oldp2)>0 and (oldp2[0] not in newrepeats):
		if commonoptions['outlog'] <= myheader.M_WARNING: print ('Warning!!! not in ', oldp2, newrepeats)

	return newrepeats

def myGMM(lendict, MinSup=2, mparameters=None, commonoptions=None):
	isillumina = False;
	if (not commonoptions==None) and commonoptions.has_key('SeqTech') and commonoptions['SeqTech']=="Illumina":
		isillumina = True;

	truecounts = None;
	if (not commonoptions==None) and commonoptions.has_key('truecounts'):
		truecounts = commonoptions['truecounts']

	minreads = int(MinSup);
	if minreads<2: minreads = 2

	ldkeys = lendict.keys(); ldkeys.sort();
	allocr = 'allocr:';
	for ldk in ldkeys:
		allocr += ('%d:%d, ' % (ldk, lendict[ldk]))
	logging.info(allocr)
	if commonoptions['outlog'] <= myheader.M_INFO: print (allocr, minreads, MinSup)

	minrepcount = 5;
	for ldk in ldkeys:
		if ldk<minrepcount: del lendict[ldk]
	ldkeys = lendict.keys(); ldkeys.sort();

	ldkeys = lendict.keys(); ldkeys.sort();
	while len(ldkeys)>1:
		lastk = ldkeys[-1]; secondlk=ldkeys[-2];
		curw = getWindowForCounts(secondlk)
		if lendict[lastk]<minreads and lastk-secondlk>curw*10:
			del lendict[lastk]
		else: break;
		ldkeys = lendict.keys(); ldkeys.sort();
	while len(ldkeys)>2:
		firstk = ldkeys[0]; secondk = ldkeys[1];
		curw = getWindowForCounts(secondk)
		if lendict[firstk]<minreads and secondk-firstk>curw*10:
			del lendict[firstk]
		else: break;
		ldkeys = lendict.keys(); ldkeys.sort();

	#for special case;
	if len(ldkeys)<1: return [[0], allocr]
	elif len(ldkeys)<2: peak2 = [ldkeys[0]];
	elif len(ldkeys)==2 or len(ldkeys)==3 or isillumina:
		maxk = ldkeys[0]
		if maxk==0: maxk = ldkeys[1]
		for ldk in ldkeys:
			if lendict[ldk]>=lendict[maxk] and ldk>0: maxk = ldk
		peak2 = [maxk]
	if len(ldkeys)<4 or isillumina:
		peak2 = checkSmallSupport(peak2, lendict, minreads)
		if len(peak2)==0: peak2 = [0]

		if isillumina:
			peak2 = getNewRepeatForIllumina(lendict, minreads, commonoptions, peak2)
		peak2.sort()
		return [peak2, allocr[:-1]];

	total_point = 0;
	mkeys = lendict.keys(); mkeys.sort();
	for mk in mkeys:
		curnum = lendict[mk] - minreads;
		if curnum<1 or mk<minrepcount: 
			continue;
		total_point += curnum
	lowcov = False;
	if total_point<50:
		total_point = 0;
		lowcov = True;
		for mk in mkeys:
			curnum = lendict[mk]
			if curnum<1 or mk<minrepcount:
				continue;
			total_point += curnum

	X = np.zeros((total_point,1))
	xi = 0;
	for mk in mkeys:
		if not lowcov:
			curnum = lendict[mk] - minreads;
		else:
			curnum = lendict[mk]
		#if curnum<1 or mk<minrepcount:
		if lendict[mk]<2 or mk<minrepcount: 
			continue;
		for j in range(curnum):
			X[xi][0] = mk; xi += 1;
	#f len(X)<200: print total_point, lowcov, len(X), X

	default_n_components = [4] #[4,3,2,5,6,4];
	#for nc in range(2, 7):
	for nc in range(3, 7):
		if nc>=total_point: break;
		if nc==4: continue
		default_n_components.append(nc)
	default_n_components.append(4)
	gmm_times = 20;
	small_covars_threhold = 0.01

	for cur_n_component_ind in range(len(default_n_components)):
		cur_n_component = default_n_components[cur_n_component_ind]
		
		atedge = True;
		for run_time in range(gmm_times):
			N = np.arange(1, (cur_n_component+1))
			models = [None for i in range(len(N))]

			for i in range(len(N)):
				models[i] = GM(N[i]).fit(X) 

			# compute the AIC and the BIC
			AIC = [m.aic(X) for m in models]
			calbic = False;
			if calbic:
				BIC = [m.bic(X) for m in models]

			mbest = models[np.argmin(AIC)]

			if commonoptions['outlog'] <= myheader.M_DEBUG:
				print ('aic', np.argmin(AIC)),
				for aic in range(len(AIC)):
					if aic==np.argmin(AIC):
						print ('<%.3f>' % (AIC[aic])),
					else:	print (' %.3f ' % (AIC[aic])),
				print ('')
				if calbic:
					print ('bic', np.argmin(BIC)),
					for bic in range(len(BIC)):
						if bic==np.argmin(BIC):
							print ('<%.3f>' % (BIC[bic])),
						else: print (' %.3f ' % (BIC[bic])),
					print ('')
			
			if 0<np.argmin(AIC)<len(AIC)-1: 
				atedge = False;
				break;
			
		has_too_small_std = False;
		for i in range(len(mbest.means_)):
			if mbest.covariances_[i,0][0]<small_covars_threhold:
				has_too_small_std = True;
		#if (not has_too_small_std) and (not atedge): break;
		if (not atedge): break;
		elif cur_n_component_ind==len(default_n_components)-1:
			if commonoptions['outlog'] <= myheader.M_WARNING:
				print ('Warning!!!! could not find optimized model')
				logging.info('Warning!!!! could not find optimized model')
			
	#print mbest.covariances_
	mean_covars = []; 
	for i in range(len(mbest.means_)):
		curk = int(mbest.means_[i,0]+0.5) #75)
		if lendict.has_key(curk):
			if commonoptions['outlog'] <= myheader.M_DEBUG: print ('>>', i, ('%9.3fm' % (mbest.means_[i,0])), ('%6d' % (lendict[curk])), ('%20.9fst' % (mbest.covariances_[i,0][0])))
			mean_covars.append([curk, mbest.means_[i,0], lendict[curk], mbest.covariances_[i,0][0]])
		else:
			closedif = sys.maxint; closekey = -1;
			for mk in mkeys:
				if closedif>abs(mk-curk):
					closedif=abs(mk-curk)
					closekey = mk
			if commonoptions['outlog'] <= myheader.M_DEBUG: print ('>>', i, ('%9.3fm' % (mbest.means_[i,0])), ('%6d' % (lendict[closekey])), ('%20.9fst' % (mbest.covariances_[i,0][0])), closekey)
			mean_covars.append([closekey, mbest.means_[i,0], lendict[closekey], mbest.covariances_[i,0][0]])

	fixed_boundarywidth = 50; close_ratio_threhold = 0.8
	
	remove_larger_covar_smaller_means = []
	for i in range(len(mean_covars)):
		mean_covars[i].append(getNeighbors_reads(lendict, mean_covars[i][0], minreads)[0])
	for i in range(len(mean_covars)):
		should_remove = False;
		for j in range(len(mean_covars)):
			if i==j: continue;
			if mean_covars[j][3]<small_covars_threhold: continue;

			if mean_covars[i][3]<=mean_covars[j][3] and mean_covars[i][1]<mean_covars[j][1]: pass
			elif mean_covars[i][3]>mean_covars[j][3] and mean_covars[i][1]<mean_covars[j][1]:
				meandif = (mean_covars[j][1]-mean_covars[i][1])*3
				#if meandif>10: meandif = 10
				if meandif>20: meandif = 20
				if mean_covars[i][3]>mean_covars[j][3]*meandif: should_remove = True
				else:
					cur_window_i = getWindowForCounts(mean_covars[i][0]) #int(mean_covars[i][0]/200.0+0.75)
					cur_window_j = getWindowForCounts(mean_covars[j][0]) #int(mean_covars[j][0]/200.0+0.75)
					if commonoptions['outlog'] <= myheader.M_INFO: print (should_remove, mean_covars[i][0], mean_covars[j][0], cur_window_i, cur_window_j, mean_covars[i][0], mean_covars[j][0], mean_covars[i][4], mean_covars[j][4], mean_covars[j][4]*close_ratio_threhold, abs(cur_window_i-cur_window_j)==1, mean_covars[i][4]<mean_covars[j][4]*close_ratio_threhold)
					
					if abs(cur_window_i-cur_window_j)==1 and abs(mean_covars[i][0]-mean_covars[j][0])<fixed_boundarywidth:
						newneighbors = getNeighbors_reads_fixed(lendict, mean_covars[i][0], cur_window_j)
						if newneighbors[0]<mean_covars[j][4]*close_ratio_threhold: should_remove = True;
					elif mean_covars[i][4]<mean_covars[j][4]*close_ratio_threhold:
						should_remove = True;
					if commonoptions['outlog'] <= myheader.M_INFO: print (should_remove)
		if not should_remove:
			remove_larger_covar_smaller_means.append(mean_covars[i])

	furtherremove = []
	for i in range(len(remove_larger_covar_smaller_means)):
		should_remove = False;
		for j in range(len(remove_larger_covar_smaller_means)):
			if i==j: continue;
			if abs(remove_larger_covar_smaller_means[i][1]-remove_larger_covar_smaller_means[j][1])<1.5 and (remove_larger_covar_smaller_means[i][3]<small_covars_threhold or remove_larger_covar_smaller_means[j][3]<small_covars_threhold):
				if remove_larger_covar_smaller_means[i][3]<small_covars_threhold and (not remove_larger_covar_smaller_means[j][3]<small_covars_threhold):
					should_remove = True;
				elif (not remove_larger_covar_smaller_means[i][3]<small_covars_threhold) and (remove_larger_covar_smaller_means[j][3]<small_covars_threhold): pass
				else:
					cur_window_i = getWindowForCounts(remove_larger_covar_smaller_means[i][0])
					cur_window_j = getWindowForCounts(remove_larger_covar_smaller_means[j][0])
					if abs(cur_window_i-cur_window_j)==1 and abs(remove_larger_covar_smaller_means[i][0]-remove_larger_covar_smaller_means[j][0])<fixed_boundarywidth:
						if cur_window_i<cur_window_j:
							newneighbors = getNeighbors_reads_fixed(lendict, remove_larger_covar_smaller_means[i][0], cur_window_j)
							if newneighbors[0]<remove_larger_covar_smaller_means[j][4]: should_remove = True;
						else:
							newneighbors = getNeighbors_reads_fixed(lendict, remove_larger_covar_smaller_means[j][0], cur_window_i)
							if newneighbors[0]>remove_larger_covar_smaller_means[i][4]: should_remove = True;
					elif remove_larger_covar_smaller_means[i][4]<remove_larger_covar_smaller_means[j][4]: should_remove = True;
		if not should_remove:
			furtherremove.append(remove_larger_covar_smaller_means[i]);

	remove_larger_covar_smaller_means, furtherremove = furtherremove, remove_larger_covar_smaller_means

	if commonoptions['outlog'] <= myheader.M_DEBUG:
		print ('keep')
		for i in range(len(remove_larger_covar_smaller_means)):
			print ('>>', i, ('%9.3fm' % (remove_larger_covar_smaller_means[i][1])), ('%6d' % (remove_larger_covar_smaller_means[i][2])), ('%20.9fst' % (remove_larger_covar_smaller_means[i][3])), mean_covars[i][4])

	peak2 = []; max_p1 = 0; max_p1_reads = 0;
	for i in range(len(remove_larger_covar_smaller_means)):
		cur_window_i = getWindowForCounts(remove_larger_covar_smaller_means[i][0]) #int(remove_larger_covar_smaller_means[i][0]/200.0+0.75)
		cur_window_max = getWindowForCounts(max_p1) #int(max_p1/200.0+0.75)
		if commonoptions['outlog'] <= myheader.M_INFO: print (max_p1, max_p1_reads, cur_window_i, cur_window_max, remove_larger_covar_smaller_means[i][0],max_p1, '<', abs(cur_window_i-cur_window_max)==1, '>', remove_larger_covar_smaller_means[i][4]>max_p1_reads)
		if abs(cur_window_i-cur_window_max)==1 and abs(remove_larger_covar_smaller_means[i][0]-max_p1)<fixed_boundarywidth:
			newneighbors = getNeighbors_reads_fixed(lendict, remove_larger_covar_smaller_means[i][0], cur_window_max) 
			if newneighbors[0] > max_p1_reads:
				max_p1 = remove_larger_covar_smaller_means[i][0]
				max_p1_reads = remove_larger_covar_smaller_means[i][4] 
		elif remove_larger_covar_smaller_means[i][4]>max_p1_reads:
			max_p1 = remove_larger_covar_smaller_means[i][0]
			max_p1_reads = remove_larger_covar_smaller_means[i][4]

	peak2.append(max_p1)
	secondpeak = {}
	for i in range(len(remove_larger_covar_smaller_means)):
		if remove_larger_covar_smaller_means[i][0]==max_p1: continue;
		if selectFromTwoX(max_p1, remove_larger_covar_smaller_means[i][0], lendict, minreads)==remove_larger_covar_smaller_means[i][0]:
			secondpeak[remove_larger_covar_smaller_means[i][0]] = []
		else:
			if remove_larger_covar_smaller_means[i][0]<max_p1 and remove_larger_covar_smaller_means[i][4]>max_p1_reads*close_ratio_threhold:
				secondpeak[remove_larger_covar_smaller_means[i][0]] = []

	if commonoptions['outlog'] <= myheader.M_INFO: print ('peak2', peak2, secondpeak)
	secondpeakkeys = secondpeak.keys();
	for spk in secondpeakkeys:
		for spkj in secondpeakkeys:	
			if spkj in [max_p1, spk]: continue;
			if selectFromTwoX(spk, spkj, lendict, minreads)==spk:
				if spk not in secondpeak[spkj]: secondpeak[spkj].append(spk)
			else:
				if spkj not in secondpeak[spk]: secondpeak[spk].append(spkj)

	for spk in secondpeakkeys:
		if  len(secondpeak[spk])==0: 
			peak2.append(spk)

	if commonoptions['outlog'] <= myheader.M_INFO: print( 'peak2', peak2, secondpeak)
	if len(secondpeakkeys)>0 and len(peak2)<2:
		max_p2 = 0; max_p2_reads = 0;
		for i in range(len(remove_larger_covar_smaller_means)):
			if secondpeak.has_key(remove_larger_covar_smaller_means[i][0]):
				if remove_larger_covar_smaller_means[i][4]>max_p2_reads: 
					max_p2_reads = remove_larger_covar_smaller_means[i][4]
					max_p2 = remove_larger_covar_smaller_means[i][0]

		peak2.append(max_p2)
	if commonoptions['outlog'] <= myheader.M_DEBUG: print ('peak2', peak2)
	peak2 = checkSmallSupport(peak2, lendict, minreads)
	if commonoptions['outlog'] <= myheader.M_DEBUG: print ('peak2', peak2)
	
	if len(peak2)==2:
		if peak2[0]==0 and (not peak2[1]==0): peak2[0] = peak2[1]
		if peak2[1]==0 and (not peak2[0]==0): peak2[1] = peak2[0]
	if len(peak2)==1 or (len(peak2)>1 and abs(peak2[0]-peak2[1])>2):
		curlen = len(peak2)
		if curlen>2: peak2 = 2;
		for i in range(curlen):
			cur_window_max = getWindowForCounts(peak2[i])
			for npdif in range(1, cur_window_max+1):
				newpeak = peak2[i]+npdif
				if lendict.has_key(newpeak) and lendict[newpeak]>1.5*lendict[peak2[i]] and lendict[newpeak]>5*minreads:
					peak2[i] = newpeak
				newpeak = peak2[i]-npdif
				if lendict.has_key(newpeak) and lendict[newpeak]>1.5*lendict[peak2[i]] and lendict[newpeak]>5*minreads:
					peak2[i] = newpeak
			#for newpeak in range(peak2[i]-cur_window_max, peak2[i]+cur_window_max+1):
			#	if lendict.has_key(newpeak) and lendict[newpeak]>2*lendict[peak2[i]]:
			#		peak2[i] = newpeak
		if len(peak2)==1:
			peak2.append(peak2[0]);
	if len(peak2)==0: peak2 = [0,0]

	peak2.sort()

	if not truecounts==None:
		if abs(truecounts[0]-peak2[0])>5 or abs(truecounts[1]-peak2[1])>5:
			if commonoptions['outlog'] <= myheader.M_INFO: print ('Big dif, please check', peak2, truecounts)

	return [peak2, allocr[:-1]]

def get2Peaks(lengd, MinSup, commonoptions=None):
	lendict = {};
	for l in lengd:
		l = int(l+0.5)
		if not lendict.has_key(l): lendict[l] = 0;
		lendict[l] += 1;

	return myGMM(lendict, MinSup, commonoptions=commonoptions)


if __name__=='__main__':
	if False:
		i=2; print (i);
		lendict = {14:24, 15:1217, 16:82, 17:2, 18:1, 19:1, 22:1, 25:1, 26:1, 31:1, 33:2, 34:1, 35:1, 50:1, 53:2, 55:1, 57:1, 60:1, 61:2, 62:1, 64:3, 65:5, 66:15, 67:27, 68:67, 69:58, 70:44, 71:21, 72:10, 73:5, 74:1, 75:1, 76:1}
		print (myGMM(lendict)[0])

		i=3; print (i);
		lendict = {20:2,43:1,44:1,45:2,46:6,47:9,48:19,49:42,50:68,51:136,52:308,53:609,54:833,55:803,56:774,57:759,58:545,59:173,60:21,61:2,62:1,64:4,65:1}
		print (myGMM(lendict)[0])

		i=4; print (i);
		lendict = {80:2,81:3,82:4,83:22,84:26,85:69,86:150,87:290,88:540,89:781,90:841,91:727,92:710,93:553,94:290,95:84,96:23,97:2,98:1,101:1}
		print (myGMM(lendict)[0])

		i=5; print (i);
		lendict = {0:85, 1:7, 2:24, 3:24, 4:30, 5:38, 6:44, 7:62, 8:52, 9:50, 10:59, 11:62, 12:61, 13:96, 14:110, 15:230, 16:461, 17:664, 18:162, 19:51, 20:4, 28:1, 30:2, 50:1, 154:1, 155:1, 215:1, 223:1, 244:1, 246:1, 329:1, 396:1, 401:1, 430:1, 463:1, 469:1, 527:1, 596:1, 600:1, 605:1, 608:1, 610:1, 613:1, 626:1, 646:1, 650:1, 655:1, 656:1, 657:1, 662:1, 664:2, 669:1, 670:1, 683:1, 685:1, 689:2, 690:1, 693:1, 695:1, 700:1, 701:1, 703:1, 707:2, 709:1, 710:1, 711:1, 715:1, 716:1, 721:2, 722:1, 723:2, 725:1, 726:2, 728:2, 731:1, 733:1, 734:1, 735:3, 740:1, 741:2, 742:1, 743:1, 744:2, 745:1, 746:2, 747:1, 748:2, 749:2, 750:2, 752:2, 753:1, 754:4, 758:2, 759:1, 760:1, 761:1, 762:3, 763:1, 764:2, 765:3, 766:4, 767:1, 768:5, 769:3, 770:3, 771:4, 772:3, 773:3, 774:6, 775:4, 776:1, 777:4, 778:3, 779:5, 780:2, 781:3, 782:8, 783:5, 784:3, 785:6, 786:4, 787:1, 788:7, 789:6, 790:2, 791:7, 792:5, 793:4, 794:10, 795:9, 796:7, 797:3, 798:5, 799:5, 800:6, 801:8, 802:6, 803:7, 804:4, 805:12, 806:8, 807:6, 808:7, 809:12, 810:10, 811:12, 812:8, 813:8, 814:10, 815:4, 816:3, 817:7, 818:12, 819:13, 820:5, 821:7, 822:12, 823:6, 824:7, 825:3, 826:6, 827:5, 828:5, 829:2, 830:5, 831:2, 832:2, 833:6, 834:2, 835:3, 836:2, 837:5, 838:3, 839:3, 840:4, 841:1, 843:2, 844:3, 845:3, 847:2, 848:2, 849:3, 850:2, 852:1, 853:2, 854:1, 855:1, 856:1, 857:1, 858:1, 860:1, 864:1, 865:1, 866:2, 868:1, 869:2, 871:2, 872:1, 875:1, 878:1, 883:2, 890:1, 893:2, 895:1, 900:1, 905:1, 908:1, 912:1, 917:1, 923:1, 952:1, 959:1, 976:1, 981:1, 1007:1, 1015:1, 1124:1}
		print (myGMM(lendict)[0])
	
		i=6; print (i);
		lendict = {0:10, 3:2, 4:4, 5:4, 6:12, 7:8, 8:9, 9:6, 10:9, 11:14, 12:8, 13:14, 14:36, 15:56, 16:161, 17:248, 18:74, 19:11, 20:1, 21:1, 26:1, 420:1, 606:1, 668:1, 680:1, 764:2, 767:1, 768:1, 775:1, 784:1, 785:2, 791:1, 793:1, 794:2, 795:1, 796:2, 797:1, 799:1, 801:1, 802:1, 803:5, 804:1, 805:3, 806:5, 808:4, 809:2, 810:4, 811:3, 812:2, 813:1, 814:6, 815:3, 816:5, 817:2, 818:3, 819:8, 820:6, 821:7, 822:7, 823:5, 824:8, 825:6, 826:8, 827:5, 828:3, 829:6, 830:6, 831:1, 832:8, 833:6, 834:6, 835:7, 837:10, 838:4, 839:3, 840:3, 841:1, 842:2, 843:2, 844:3, 845:3, 846:3, 847:3, 848:4, 849:3, 850:2, 851:4, 852:2, 853:1, 856:1, 858:1, 859:2, 860:2, 862:1, 863:1, 864:2, 867:1, 869:1, 872:1, 874:1, 875:1, 877:1, 885:1, 890:1, 894:1, 895:1, 896:1, 897:1, 910:1, 915:1, 921:1, 922:1, 924:1, 996:1, 1023:1, 1074:1, 1121:1}
		print (myGMM(lendict)[0])
	
		i=7; print (i);
		lendict = {0:56, 1:10, 2:39, 3:44, 4:33, 5:57, 6:58, 7:60, 8:60, 9:71, 10:67, 11:71, 12:78, 13:113, 14:195, 15:452, 16:976, 17:1170, 18:423, 19:78, 20:21, 21:4, 22:1, 24:1, 33:1, 57:1, 63:1, 95:1, 96:1, 97:1, 144:1, 157:1, 161:1, 189:1, 197:1, 211:1, 215:1, 219:1, 221:1, 228:1, 229:1, 240:1, 242:1, 243:1, 245:1, 249:1, 260:1, 271:1, 278:1, 287:1, 292:1, 296:1, 301:1, 309:1, 331:1, 339:1, 343:1, 344:1, 352:1, 355:1, 357:1, 358:1, 369:2, 370:1, 372:1, 374:2, 375:1, 378:1, 379:1, 381:1, 384:1, 385:1, 387:2, 388:1, 389:1, 390:1, 392:1, 393:1, 395:2, 396:2, 397:1, 398:1, 399:1, 400:1, 401:1, 402:2, 403:1, 405:1, 406:3, 407:1, 408:2, 409:1, 410:4, 411:2, 412:3, 413:2, 414:1, 415:1, 416:3, 418:3, 419:2, 420:2, 421:3, 422:2, 423:2, 424:2, 425:3, 426:1, 427:3, 428:4, 429:6, 430:2, 431:2, 432:4, 433:6, 434:6, 435:5, 436:1, 437:7, 438:6, 439:3, 440:8, 441:4, 442:7, 443:11, 444:4, 445:9, 446:6, 447:5, 448:6, 449:4, 450:3, 451:9, 452:12, 453:6, 454:7, 455:6, 456:7, 457:11, 458:10, 459:6, 460:11, 461:5, 462:8, 463:12, 464:4, 465:6, 466:6, 467:15, 468:12, 469:15, 470:11, 471:14, 472:17, 473:11, 474:15, 475:15, 476:20, 477:9, 478:13, 479:13, 480:21, 481:14, 482:16, 483:13, 484:18, 485:23, 486:14, 487:8, 488:12, 489:13, 490:8, 491:10, 492:4, 493:7, 494:4, 495:4, 496:8, 497:3, 498:8, 499:3, 500:2, 501:5, 502:6, 503:3, 504:3, 505:3, 506:3, 508:2, 510:2, 511:3, 513:1, 514:2, 518:1, 524:1, 525:2, 528:1, 530:1, 534:1, 543:1, 544:1, 549:1, 553:1, 554:1, 556:2, 563:1, 564:1, 567:1, 577:1, 583:1, 595:1, 603:1, 616:1, 626:1, 628:1, 643:1, 646:1, 680:1, 732:1, 826:1}
		print (myGMM(lendict)[0])

		i=8; print (i);
		lendict = {431:2, 432:4, 433:6, 434:6, 435:5, 436:1, 437:7, 438:6, 439:3, 440:8, 441:4, 442:7, 443:11, 444:4, 445:9, 446:6, 447:5, 448:6, 449:4, 450:3, 451:9, 452:12, 453:6, 454:7, 455:6, 456:7, 457:11, 458:10, 459:6, 460:11, 461:5, 462:8, 463:12, 464:4, 465:6, 466:6, 467:15, 468:12, 469:15, 470:11, 471:14, 472:17, 473:11, 474:15, 475:15, 476:20, 477:9, 478:13, 479:13, 480:21, 481:14, 482:16, 483:13, 484:18, 485:23, 486:14, 487:8, 488:12, 489:13, 490:8, 491:10, 492:4, 493:7, 494:4, 495:4, 496:8, 497:3, 498:8, 499:3, 500:2, 501:5, 502:6, 503:3, 504:3, 505:3, 506:3, 508:2, 510:2, 511:3}
		print myGMM(lendict)[0]

		i=9; print (i);
		lendict = {0:64, 3:40, 4:68, 5:56, 6:130, 7:184, 8:255, 9:177, 10:313, 11:371, 12:755, 13:1540, 14:1608, 15:642, 16:200, 17:53, 18:16, 19:11, 20:5, 21:1, 22:3, 23:1, 24:1}
		print myGMM(lendict)[0]

		i=10; print (i);
		lendict = {0:700, 1:2, 2:3, 3:13, 4:15, 5:14, 6:18, 7:13, 8:23, 9:29, 10:21, 11:49, 12:52, 13:59, 14:89, 15:106, 16:132, 17:200, 18:296, 19:351, 20:283, 21:158, 22:74, 23:28, 24:9, 25:11, 26:16, 27:9, 28:5, 29:5, 30:7, 31:7, 32:5, 33:8, 34:9, 35:5, 36:6, 37:5, 38:11, 39:11, 40:3, 41:6, 42:10, 43:12, 44:6, 45:21, 46:23, 47:17, 48:22, 49:39, 50:42, 51:87, 52:141, 53:197, 54:259, 55:427, 56:568, 57:809, 58:1133, 59:1485, 60:1930, 61:2184, 62:2211, 63:1808, 64:1266, 65:750, 66:429, 67:243, 68:143, 69:112, 70:80, 71:59, 72:43, 73:29, 74:27, 75:14, 76:15, 77:10, 78:15, 79:12, 80:6, 81:4, 82:4, 83:5, 84:1, 85:3, 86:2, 87:2, 88:3, 89:3, 90:1, 91:1, 92:2, 93:1, 96:1, 97:1, 99:1, 101:1, 104:1, 106:2, 109:1, 130:1, 148:1, 153:1, 188:1}
		print (myGMM(lendict)[0])

		i=11; print (i);
		lendict = {50:1, 51:1, 52:2, 53:1, 54:2, 55:2, 56:4, 57:5, 58:1}
		print (myGMM(lendict)[0])

	else:
		lendict = {38:7, 39:14, 40:33, 41:33, 42:68, 43:159, 44:445, 45:856, 46:790, 47:326, 48:91, 49:46, 50:99, 51:238, 52:508, 53:671, 54:515, 55:179, 56:31, 57:6, 59:2}
		
		lendict = {16:1, 37:3, 38:1, 39:5, 40:12, 41:32, 42:30, 43:66, 44:161, 45:482, 46:908, 47:746, 48:282, 49:90, 50:46, 51:79, 52:226, 53:544, 54:717, 55:474, 56:174, 57:36, 58:3}

		lendict = {18:2, 25:2, 26:2, 27:2, 28:6, 29:7, 30:4, 31:30, 32:54, 33:91, 34:72, 35:106, 36:150, 37:510, 38:1018, 39:838, 40:278, 41:49, 42:4, 56:2, 57:1, 58:2, 59:2, 60:3, 61:5, 62:16, 63:29, 64:47, 65:67, 66:65, 67:116, 68:260, 69:446, 70:376, 71:255, 72:63, 73:16, 75:1}
		
		lendict = {0:362, 1:3, 2:5, 3:12, 4:13, 5:11, 6:20, 7:13, 8:32, 9:54, 10:56, 11:51, 12:81, 13:61, 14:98, 15:92, 16:99, 17:106, 18:144, 19:173, 20:208, 21:219, 22:334, 23:751, 24:1690, 25:2847, 26:3590, 27:2586, 28:830, 29:429, 30:463, 31:780, 32:1396, 33:1947, 34:1954, 35:1214, 36:436, 37:146, 38:59, 39:31, 40:20, 41:24, 42:6, 43:13, 44:4, 45:3, 46:2, 47:1, 51:2, 52:1, 53:1, 58:1, 72:1}

		lendict = {11:40, 12:36, 13:37, 14:37, 15:74, 16:41, 17:35, 18:32, 19:27, 20:31, 21:17, 22:23, 23:24, 24:22, 25:26, 26:25, 27:18, 28:12, 29:17, 30:12, 31:8, 55:3, 56:3, 57:4, 58:5, 59:6, 60:14, 61:16, 62:13, 63:15, 64:8, 65:5, 66:2}	
		lendict = {0:35, 2:3, 3:14, 4:10, 5:13, 6:21, 7:39, 8:49, 9:67, 10:88, 11:146, 12:428, 13:1309, 14:2805, 15:1006, 16:273, 17:107, 18:38, 19:16, 20:15, 21:4, 22:1, 23:2, 24:3, 26:2, 31:1}
		lendict = {6:2, 7:3, 8:2, 9:5, 10:8, 11:13, 12:4, 13:7, 14:11, 15:8, 16:1, 19:1}

		commonoptions = {}
		commonoptions['outlog'] = myheader.M_DEBUG

		lendict = {15:1, 16:1, 17:1, 22:1, 27:1, 31:1, 34:1, 35:1, 38:1, 41:1, 47:1, 48:1, 49:1, 51:2, 52:2, 53:2, 54:3, 55:2, 56:2, 58:3, 59:2, 60:3, 61:4, 62:1, 63:2, 64:2, 65:3, 66:1, 68:2, 69:1, 70:1, 71:1, 72:2, 73:2, 75:4, 76:1, 78:2, 79:4, 80:1, 81:1, 82:2, 83:1, 84:1, 87:2, 89:1, 90:1, 94:2, 95:1, 96:1, 97:1, 100:1, 101:1, 103:1, 107:1, 1318:1}
		lendict = {0:8, 1:2, 2:14, 3:23, 4:19, 5:29, 6:31, 7:40, 8:52, 9:41, 10:8, 16:2, 19:1, 24:1, 26:1}

		lendict = {11:3, 12:32, 13:469, 14:3967, 15:8685, 16:5655, 17:2141, 18:931, 19:447, 20:285, 21:174, 22:109, 23:84, 24:69, 25:50, 26:29, 27:43, 28:19, 29:18, 30:14, 31:11, 32:14, 33:8, 34:12, 35:5, 36:5, 37:10, 38:6, 39:3, 40:7, 41:1, 42:4, 43:4, 44:1, 45:4, 46:4, 47:1, 48:2, 49:1, 50:1, 52:1, 53:3, 56:2, 57:1, 58:2, 61:1, 62:1, 64:2, 65:1, 66:3, 68:1, 69:1, 70:2, 71:1, 72:4, 73:3, 74:2, 75:3, 76:4, 77:2, 79:3, 80:2, 81:5, 82:15, 83:12, 84:20, 85:30, 86:48, 87:59, 88:103, 89:124, 90:147, 91:131, 92:148, 93:136, 94:120, 95:104, 96:82, 97:57, 98:55, 99:42, 100:27, 101:43, 102:24, 103:15, 104:15, 105:11, 106:9, 107:9, 108:8, 109:6, 110:11, 111:8, 112:7, 113:5, 114:4, 115:3, 116:4, 117:5, 118:3, 119:5, 120:6, 121:1, 122:1, 123:3, 124:1, 125:3, 126:1, 128:1, 129:1, 130:1, 131:1, 134:1, 135:1, 136:2, 137:1, 140:1, 143:1, 144:1, 147:1, 148:1, 149:1, 152:1, 155:1, 162:1, 169:1, 171:2, 185:1, 204:1, 235:1, 240:1, 316:1, 376:1, 386:1, 429:1, 444:1, 520:1, 532:1, 557:1, 560:1, 571:1, 572:1, 573:1, 574:1, 578:1, 580:1, 592:1, 593:1, 604:1, 612:1, 620:1, 622:1, 626:1, 632:1, 640:1, 643:1, 661:1, 679:1, 694:1, 727:1, 736:1, 758:1, 772:1, 776:1, 786:2, 797:1, 806:1, 835:1, 866:1, 868:1, 906:1, 1902:1}

		print (myGMM(lendict, commonoptions=commonoptions)[0])
