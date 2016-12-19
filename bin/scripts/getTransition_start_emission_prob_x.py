
import os;
import sys;

import string

import numpy as np;

import getTransition_start_emission_prob_without0
import printHMMmatrix

import getTransition_start_emission_prob_2
import getTransition_start_emission_prob_3
import getTransition_start_emission_prob_4
import getTransition_start_emission_prob_5
import getTransition_start_emission_prob_6

import getTransition_start_emission_prob_x_without0

from myheader import *

def getTransition_start_emission_prob_x(repPat, forprint=False):
	if len(repPat)<1: return None

	avgsub = 0.0005
	avgsub = hmm_random_rep_transit/len(repPat)
	repPat = string.strip(repPat);
	
	typeOfRepEle = ['', 'I', 'D'];
	repEle = [];
	for rp_ind in range(len(repPat)):
		repEle.append(''.join(['r', str(rp_ind+1)]));
	states = ['N'];
	for typRE in typeOfRepEle:
		for rp in repEle:
			states.append(''.join([typRE, rp]));

	trainsmat = np.full((len(states), len(states)), 1e-9);	
	#for N to N
	trainsmat[0][0] = 0.96;
	#for N to rep;
	if not len(repPat)<2:
		trainsmat[0][1] = 0.02;
	else: trainsmat[0][1] = 0.04;
	if not len(repPat)<2:
		trainsmat[0][1+len(repEle)*2] = 0.02;
	#for rep to N;
	trainsmat[len(repEle)][0] = 0.02;
	trainsmat[len(repEle)*2][0] = 0.02;
	trainsmat[len(repEle)*3][0] = 0.02;
	if not len(repPat)<2:
		trainsmat[len(repEle)*3-1][0] = 0.02;
	#avgsub
	for i in range(1, len(states)):
		for j in range(len(repEle)):
			trainsmat[i][j+1] = avgsub
	#for insertion
	add_index = len(repEle)+1;
	for typ_ind in range(len(typeOfRepEle)):
		for j in range(len(repEle)):
			if typ_ind<len(typeOfRepEle)-1:
				jind = j
			else:
				jind = j+1;
				if jind > len(repEle)-1:
					jind = 0;
			trainsmat[len(repEle)*typ_ind+j+1][jind+add_index] = 0.11
	#for deletion
	add_index = len(repEle)*2+1;
	for typ_ind in range(len(typeOfRepEle)):
		for j in range(len(repEle)):
			if len(repPat)<2: continue;
			if typ_ind<len(typeOfRepEle)-1:
				jind = j+1;
				if jind > len(repEle)-1:
					jind = 0;
			else:
				jind = j+2
				if jind > len(repEle)-1:
					jind -= len(repEle)
			trainsmat[len(repEle)*typ_ind+j+1][jind+add_index] = 0.02
	#for between rep
	add_index = 1;
	for typ_ind in range(len(typeOfRepEle)):
		for j in range(len(repEle)):
			if typ_ind<len(typeOfRepEle)-1:
				jind = j+1;
				if jind > len(repEle)-1:
					jind = 0;
			else:
				jind = j+2
				if jind > len(repEle)-1:
					jind -= len(repEle)
			restprob = 1;
			for jstat in range(len(states)):
				if jstat==jind+add_index: pass
				else: restprob -= trainsmat[len(repEle)*typ_ind+j+1][jstat];
			trainsmat[len(repEle)*typ_ind+j+1][jind+add_index] = restprob
	
	startprob = []
	for i in range(len(states)):
		startprob.append(1e-9)
	startprob[0] = 0.96;
	if len(repPat)<2:
		startprob[1] = 0.04
	else:
		startprob[1] = 0.02
		startprob[1+len(repEle)*2] = 0.02
	startprob = np.array(startprob)

	emisionmat = np.full((len(repEle)*len(typeOfRepEle)+1, 4), 0.005);
	randrow = [0]
	for j in range(len(repEle)):
		randrow.append(j+len(repEle)+1);
	if len(repPat)<2: randrow.append(len(repEle)*len(typeOfRepEle))
	#print randrow
	for rdr in randrow:
		for jcol in range(4):
			emisionmat[rdr][jcol] = 0.25;

        obs_symbols = np.array(['A', 'C', 'G', 'T'])
        for naind in range(len(repPat)):
                emind = (np.where(obs_symbols==repPat[naind]))[0][0]
                emisionmat[naind+1][emind] = 0.985
		if len(repPat)<2: continue;
                if naind<len(repPat)-1:
                        afterd = naind + 1;
                else:
                        afterd = 0;
                emind = (np.where(obs_symbols==repPat[afterd]))[0][0]
                emisionmat[naind+1+len(repPat)*2][emind] = 0.985

        if getTransition_start_emission_prob_without0.outputm or forprint:
                print 'HMMmatrix1'
                printHMMmatrix.printHMMmatrix(states, obs_symbols, trainsmat, emisionmat, startprob)

        state3class = [range(1, len(repPat)+1), range(len(repPat)+1, 2*len(repPat)+1), range(2*len(repPat)+1, 3*len(repPat)+1)]

        return [trainsmat, startprob, emisionmat, obs_symbols, np.array(states), len(states), len(obs_symbols), state3class]

def compareTwoFloat(a, b):
	if isinstance(a, float) and isinstance(b, float):
		if int(a*10000)==int(b*10000) or int(a*10000)==int(b*10000+0.5) or int(a*10000+0.5)==int(b*10000): return True;
		else: return False;
	else: return a==b;

def compareTwoNumpyArray(a1, a2):
	isSame = True;
	for i in range(a1.shape[0]):
		for j in range(a1.shape[1]):
			if compareTwoFloat(a1[i][j], a2[i][j]): pass;
			else: 
				isSame = False;
				print ('\t <%d,%d> 1=%.10f, 2=%.10f' % (i, j, a1[i][j], a2[i][j]))
	return isSame
def compareTwoNumpyArray1(a1, a2):
	isSame = True;
	for i in range(a1.shape[0]):
		if compareTwoFloat(a1[i], a2[i]): pass;
		else:
			isSame = False;
			print ('\t <%d> 1=%.10f, 2=%.10f' % (i, a1[i], a2[i]))
	return isSame
def compareTwoMat(a1, a2):
        isSame = True;
        for i in range(len(a1)):
                for j in range(len(a1[0])):
                        if compareTwoFloat(a1[i][j], a2[i][j]): pass;
                        else:
                                isSame = False;
                                print ('\t <%d,%d> 1=%.10f, 2=%.10f' % (i, j, a1[i][j], a2[i][j]))
        return isSame
def compareTwoMat1(a1, a2):
	isSame = True;
	for i in range(len(a1)):
		if compareTwoFloat(a1[i], a2[i]): pass;
		else:
			isSame = False;
			#print ('\t <%d> 1=%.10f, 2=%.10f' % (i, a1[i], a2[i]))
			print ('\t <%d> 1=%s, 2=%s' % (i, a1[i], a2[i]))
	return isSame
def CompareTwoNumpyArrays(matnames, matindex, fun, matorg, matx):
	for mid in range(len(matnames)):
		a1 = matorg[matindex[mid]];  a2 = matx[matindex[mid]]
		print matnames[mid],
		if fun(a1, a2):
			print ': Same'

def compareMat(matorg, matx):
	matnames = ['trainsmat', 'emisionmat']
	matindex = [    0,          2]
	CompareTwoNumpyArrays(matnames, matindex, compareTwoNumpyArray, matorg, matx)

	matnames = ['startprob', 'obs_symbols']
	matindex = [   1,            3,       ]
	CompareTwoNumpyArrays(matnames, matindex, compareTwoNumpyArray1, matorg, matx)

	CompareTwoNumpyArrays(['states'], [4], compareTwoMat1, matorg, matx)

	if matorg[5]==matx[5]: print 'state_num is same'
	else: print ('state_num is not same: %d, %d' % (matorg[5], matx[5]))
	if matorg[6]==matx[6]: print 'state_num is same'
	else: print ('symbols_num is not same: %d, %d' % (matorg[6], matx[6]))
	
	print 'state3class',
	if compareTwoMat(matorg[7], matx[7]):
		print ': Same'
	

if __name__=='__main__':
	wouldprint = False; #True;

	print 'len(Pattern) = 4'
	matorg4 = getTransition_start_emission_prob_4.getTransition_start_emission_prob_4('TATC', wouldprint)
	matx4   = getTransition_start_emission_prob_x('TATC', wouldprint);
	compareMat(matorg4, matx4)

	print '\nlen(Pattern) = 3'
	matorg3 = getTransition_start_emission_prob_3.getTransition_start_emission_prob_3('CAG', wouldprint)
	matx3   = getTransition_start_emission_prob_x('CAG', wouldprint);
	compareMat(matorg3, matx3)

	print '\nlen(Pattern) = 2'
	matorg2 = getTransition_start_emission_prob_2.getTransition_start_emission_prob_2('CG', wouldprint)
	matx2   = getTransition_start_emission_prob_x('CG', wouldprint);
	compareMat(matorg2, matx2)
	
	print '\nlen(Pattern) = 5'
	matorg5 = getTransition_start_emission_prob_5.getTransition_start_emission_prob_5('TATCG', wouldprint)
	matx5   = getTransition_start_emission_prob_x('TATCG', wouldprint);
	compareMat(matorg5, matx5)

	print '\nlen(Pattern) = 6'
	matorg6 = getTransition_start_emission_prob_6.getTransition_start_emission_prob_6('TATCGG', wouldprint)
	matx6   = getTransition_start_emission_prob_x('TATCGG', wouldprint);
	compareMat(matorg6, matx6)

	wouldprint = True; 
	getTransition_start_emission_prob_x_without0.getTransition_start_emission_prob_x_without0('CAG', wouldprint)
	
	matx1   = getTransition_start_emission_prob_x('G', wouldprint);

	#etTransition_start_emission_prob_without0.getTransition_start_emission_prob_without0('TATC')

