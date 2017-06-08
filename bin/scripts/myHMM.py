
import string

import numpy as np;
from hmmlearn import hmm

import logging

from myheader import *

import getTransition_start_emission_prob_x
import printHMMmatrix

def getBasePair():
        bp = {}
        bp['A'] = 'T'
        bp['C'] = 'G'
        bp['G'] = 'C'
        bp['T'] = 'A'

        return bp;
def getComplementary(bp, na):
        return bp[na.upper()]

def getComplementary3(bp, na3):
        c3 = ''
        for li in range(len(na3)):
                c3 += getComplementary(bp, na3[li]);
        return (c3[::-1]);

def getComplementaryCompPat(bp, commonOptions):
	revs = []
	for i in range(len(commonOptions['CompRep'])):
		revs.append({})
		compkeys = commonOptions['CompRep'][i].keys();
		for ck in compkeys:
			revs[-1][getComplementary(bp,ck)] = commonOptions['CompRep'][i][ck]
	commonOptions['CompRep'] = revs[::-1]
	
def produce_for_repPat(commonOptions, moreOptions):
	forw_rerv = moreOptions['forw_rerv']
	repPat = moreOptions['repPat']
	bp = getBasePair()
	if forw_rerv[0]=='-':
		if commonOptions['CompRep']=='0': 
			moreOptions['repPat'] = getComplementary3(bp, repPat)
		else:
			getComplementaryCompPat(bp, commonOptions)
		moreOptions['forw_rerv'] ='+' + moreOptions['forw_rerv'][1:]

	if not commonOptions['CompRep']=='0':
		repplist = []
		for i in range(len(commonOptions['CompRep'])):
			compkeys = commonOptions['CompRep'][i].keys(); compkeys.sort()
			repplist.append(compkeys[0])
			for ck in compkeys:
				if commonOptions['CompRep'][i][ck]>commonOptions['CompRep'][i][repplist[-1]]:
					repplist[-1] = ck
		moreOptions['repPat'] = ''.join(repplist)

	if commonOptions['outlog'] <= M_INFO:
		print 'produce_for_repPat', moreOptions['repPat'], commonOptions['CompRep']
		logging.info('produce_for_repPat ' + str(moreOptions['repPat']) +' ' + str(commonOptions['CompRep']))

def getTransition_start_emission_prob(repPat, commonOptions, forprint=False):
	return getTransition_start_emission_prob_x.getTransition_start_emission_prob_x(repPat, commonOptions, forprint);

#coded on 11 Jan 2017
def getPred(predstats, obs_seq, state3class, patlen):
	newstr = []; ststar = []; stategrouplist = [] 
	for stind_ind in range(len(predstats)):
		stind = predstats[stind_ind]
		ststar.append(str(stind))
	
		if stind==0:
			if len(newstr)==0 or stind_ind==0 or (not predstats[stind_ind-1]==0):
				newstr.append([]); stategrouplist.append(0);
			newstr[-1].append(obs_seq[stind_ind]);
		elif stind in state3class[0]:
			if len(newstr)==0 or stind_ind==0 or (predstats[stind_ind-1]==0):
				newstr.append([]); stategrouplist.append(1);
			newstr[-1].append(obs_seq[stind_ind])
		elif stind in state3class[1]: 
			pass
		elif stind in state3class[2]:
			if len(newstr)==0 or stind_ind==0 or (predstats[stind_ind-1]==0):
				newstr.append([]); stategrouplist.append(1);
			#
			if stind_ind==0:
				pre_stind = len(state3class[2])
			else:
				pre_stind = (predstats[stind_ind-1])
			cur_mod = (stind-1)%len(state3class[2])
			if pre_stind>2*len(state3class[2]):
				pre_mod = (pre_stind)%len(state3class[2])
			else:
				pre_mod = (pre_stind-1)%len(state3class[2])
			for i in range(0, len(state3class[2])):
				next_mod = pre_mod + i;
				if next_mod>=len(state3class[2]): next_mod -= len(state3class[2])
				if next_mod==cur_mod: break;
				newstr[-1].append('-')
				if i>0: 
					del_str = stind_ind-3; 
					if del_str<0: del_str = 0
					del_end = stind_ind+3
					if del_end>=len(predstats): del_end = len(predstats)-1
					if cur_M_STAT <= M_INFO: print 'mutiple deletion: ', stind_ind, pre_stind, stind, predstats[del_str:del_end], obs_seq[del_str:del_end]
			#
			newstr[-1].append(obs_seq[stind_ind])
		else:
			if cur_M_STAT <= M_WARNING: print 'Warning!!! in hmm Wrong state '+str(stind)


	#remove isolated regions
	large_start_ind = 0; large_end_ind = len(newstr)-1
	while large_start_ind<=large_end_ind:
		if stategrouplist[large_start_ind]==0: 
			large_start_ind += 1; 
			continue;
		if stategrouplist[large_end_ind]==0: 
			large_end_ind -= 1
			continue;

		rem = False;
		if len(newstr[large_start_ind])>=len(newstr[large_end_ind]):
			large_start_ind, large_end_ind, currem2 = add_end(large_start_ind, large_end_ind, newstr, predstats, patlen, stategrouplist)
			large_start_ind, large_end_ind, currem1 = add_start(large_start_ind, large_end_ind, newstr, predstats, patlen, stategrouplist)
		else:
			large_start_ind, large_end_ind, currem1 = add_start(large_start_ind, large_end_ind, newstr, predstats, patlen, stategrouplist)
			large_start_ind, large_end_ind, currem2 = add_end(large_start_ind, large_end_ind, newstr, predstats, patlen, stategrouplist)
		
		if currem1 or currem2: rem = True;

		if not rem: break

	pre0 = 0; large_rep = ['']
	if large_start_ind<=large_end_ind:
		for i in range(large_start_ind):
			pre0 += len(newstr[i])
		large_rep.extend(newstr[large_start_ind])
		large_start_ind += 1
		while large_start_ind<=large_end_ind:
			large_rep.extend(newstr[large_start_ind])
			large_start_ind += 1
	else: large_rep = ['']

	return [''.join(large_rep), ''.join(ststar), pre0]

def add_start(large_start_ind, large_end_ind, newstr, predstats, patlen, stategrouplist, repover0=0.7): #repover0=0.7): #repover0=2):
	rem = False;
	if large_start_ind+1<large_end_ind:
		cur_rep_len = len(newstr[large_start_ind])
		if not stategrouplist[large_start_ind+1]==0:
			if cur_M_STAT <= M_WARNING: print 'In hmm Wrong split f', predstats, newstr, large_start_ind, large_start_ind+1, large_end_ind
		next_0_len = len(newstr[large_start_ind+1])
		if cur_rep_len/float(next_0_len)<repover0 and cur_rep_len/float(patlen)<len_isolated_repeat:
			large_start_ind += 2;
			rem = True;
	return [large_start_ind, large_end_ind, rem]

def add_end(large_start_ind, large_end_ind, newstr, predstats, patlen, stategrouplist, repover0=0.7): #repover0=0.7): #repover0=2):
	rem = False;
	if large_start_ind<large_end_ind-1:
		cur_rep_len = len(newstr[large_end_ind])
		if not stategrouplist[large_end_ind-1]==0:
			if cur_M_STAT <= M_WARNING: print 'In hmm Wrong split b', predstats, newstr, large_start_ind, large_end_ind-1, large_end_ind
		next_0_len = len(newstr[large_end_ind-1])
		if cur_rep_len/float(next_0_len)<repover0 and cur_rep_len/float(patlen)<len_isolated_repeat:
			large_end_ind -= 2
			rem = True;
	return [large_start_ind, large_end_ind, rem]

def hmmpred(obs_seq, na3, forw_rerv, hmmoptions, commonOptions):
	obs_seq = obs_seq.replace('-', '')
	#obs_seq = obs_seq.replace('N', ''); obs_seq = obs_seq.replace('n', '');
	
	bp = getBasePair()
	len_repPat = printHMMmatrix.get_len_repPat(na3, commonOptions)

	trainsmat, startprob, emisionmat, obs_symbols, states, numStates, numSymbols, state3class, tol_info = hmmoptions
	hmmmodel = hmm.MultinomialHMM(numStates)
	hmmmodel.transmat_ = trainsmat;
	hmmmodel.startprob_ = startprob;
	hmmmodel.emissionprob_ = emisionmat;
	hmmmodel.n_features = numSymbols

	myobs = []
	for osi in range(len(obs_seq)): 
		myobs.append((np.where(obs_symbols==obs_seq[osi]))[0][0]) 

	logprob, predstats = hmmmodel.decode(np.array([myobs]).T, algorithm="viterbi")

	newstr, ststar, pre0 = getPred(predstats, obs_seq, state3class, len_repPat)
	if cur_M_STAT <= M_DEBUG: #int(len(newstr)/float(len_repPat)+0.5)<14: #False: #True: #False: #int(len(newstr)/float(len_repPat)) in [8,13]:
		print 'hmmB:', obs_seq, int(len(newstr)/float(len_repPat)+0.5)
		psstr = []
		for ps in predstats: psstr.append(str(ps))
		print 'hmmB:', ''.join(psstr)

	return [newstr, pre0, ststar]

if __name__=='__main__':
	na3 = 'CTG'; forw_rerv='+'
	obs_seq = 'ATAGGTCCCCCTGCTGCTGCTGCTGCTGCTGCTGTTGCTGCTTTTGCTGATGTCTGAAAC'
	#obs_seq = 'ATAGTCCCCCTGCTGCATCTGCTGTAGCTGCATGTAAACATAGTAAGCAATATAATAATAAAAAATATATAAGGAAAAACA'
	obs_seq = 'CTGAGTGTTGAGCTGGCTGACTAGATTGCTTGGTATCGTGCTTTTGCGCTG'

	##na3= 'CAG'
	##obs_seq = 'CTGATGTTGCTGCTGCTGATGCTGCTGCTGCGCTG'
	##print hmmpred(obs_seq, na3, forw_rerv, 0)
	##obs_seq = 'CTTGTCCGATATGTTCTCTCGCTGTGGCTTCTTGCTCTTGGGGCTAAGTGGATTCTGTATGACTGTCGAACTACATGCAAAAAGT'

	obs_seq = 'ATAGGTCCCCCTGCTGCTGCTGCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTGCTGCTGTTGCTGCTTTTGCTGATGTCTGAAAC'
	print hmmpred(obs_seq, na3, forw_rerv, getTransition_start_emission_prob(na3))


