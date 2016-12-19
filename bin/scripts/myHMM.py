
import string

import numpy as np;
from hmmlearn import hmm

#import getTransition_start_emission_prob_2
#import getTransition_start_emission_prob_3
#import getTransition_start_emission_prob_4
#import getTransition_start_emission_prob_5
#import getTransition_start_emission_prob_6
#import getTransition_start_emission_prob_without0

import getTransition_start_emission_prob_x
import getTransition_start_emission_prob_x_without0

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

def getTransition_start_emission_prob(repPat, forprint=False):
	return getTransition_start_emission_prob_x.getTransition_start_emission_prob_x(repPat, forprint);
	
	#if len(repPat)==2: return   getTransition_start_emission_prob_2.getTransition_start_emission_prob_2(repPat, forprint);
	#elif len(repPat)==3: return getTransition_start_emission_prob_3.getTransition_start_emission_prob_3(repPat, forprint);
	#elif len(repPat)==4: return getTransition_start_emission_prob_4.getTransition_start_emission_prob_4(repPat, forprint);
	#elif len(repPat)==5: return getTransition_start_emission_prob_5.getTransition_start_emission_prob_5(repPat, forprint);
	#elif len(repPat)==6: return getTransition_start_emission_prob_6.getTransition_start_emission_prob_6(repPat, forprint);
	#else: return None 

def getTransition_start_emission_prob_without0(repPat, forprint=False):
	return getTransition_start_emission_prob_x_without0.getTransition_start_emission_prob_x_without0(repPat, forprint)
	
	#return getTransition_start_emission_prob_without0.getTransition_start_emission_prob_without0(repPat, forprint)

def getPred(predstats, obs_seq, state3class, state_add=0):
        newstr = ''; ststar = ''; pre0 = 0; ispre = True;
        for stind_ind in range(len(predstats)):
                stind = predstats[stind_ind] + state_add
                ststar += str(stind)
                if stind==0:
			if ispre: pre0 += 1;
		else: ispre = False;
                #if stind in [1,2,3]: newstr += obs_seq[stind_ind] #(states[stind]);
                #if stind in [4,5,6]: pass
                #if stind in [7,8,9]:
                if stind in state3class[0]: newstr += obs_seq[stind_ind] #(states[stind]);
                if stind in state3class[1]: pass
                if stind in state3class[2]: 
			newstr += ('-'); 
			newstr += obs_seq[stind_ind]
 	return [newstr, ststar, pre0]

def hmmpred(obs_seq, na3, forw_rerv, hmmoptions, afterbefore = 1):
	bp = getBasePair()
	if forw_rerv[0]=='-': na3 = getComplementary3(bp, na3)
	#if afterbefore>0:
	#	trainsmat, startprob, emisionmat, obs_symbols, states, numStates, numSymbols = getTransition_start_emission_prob(na3)
	#else:
	#	trainsmat, startprob, emisionmat, obs_symbols, states, numStates, numSymbols = getTransition_start_emission_prob_without0(na3)
	trainsmat, startprob, emisionmat, obs_symbols, states, numStates, numSymbols, state3class = hmmoptions
	hmmmodel = hmm.MultinomialHMM(numStates)
	hmmmodel.transmat_ = trainsmat;
	hmmmodel.startprob_ = startprob;
	hmmmodel.emissionprob_ = emisionmat;
	hmmmodel.n_features = numSymbols

	myobs = []
	for osi in range(len(obs_seq)): myobs.append((np.where(obs_symbols==obs_seq[osi]))[0][0]) 

	logprob, predstats = hmmmodel.decode(np.array([myobs]).T, algorithm="viterbi")

	#print logprob
	#print predstats	
	#0     1    2    3     4    5     6      7    8     9
	#'N', 'C', 'T', 'G', 'IC', 'IT', 'IG', 'DC', 'DT', 'DG'
	if afterbefore>0:
		newstr, ststar, pre0 = getPred(predstats, obs_seq, state3class)
	else: newstr, ststar, pre0 = getPred(predstats, obs_seq, state3class, 1)

	#print numStates, numSymbols, state3class, newstr, ststar, pre0

	moutput = False;
	if moutput: # or len(obs_seq)>120:
		prestr = '';
		for i in range(pre0): prestr += ' ';	
		print na3, ': o2', obs_seq, '\n', na3, ': sp', ststar, '\n', na3, ': p=',  prestr+newstr, '\t', len(newstr)	

	return [newstr, pre0, ststar]

if __name__=='__main__':
	na3 = 'CTG'; forw_rerv='+'
	obs_seq = 'ATAGGTCCCCCTGCTGCTGCTGCTGCTGCTGCTGTTGCTGCTTTTGCTGATGTCTGAAAC'
	#obs_seq = 'ATAGTCCCCCTGCTGCATCTGCTGTAGCTGCATGTAAACATAGTAAGCAATATAATAATAAAAAATATATAAGGAAAAACA'
	obs_seq = 'CTGAGTGTTGAGCTGGCTGACTAGATTGCTTGGTATCGTGCTTTTGCGCTG'

	#na3= 'CAG'
	#obs_seq = 'CTGATGTTGCTGCTGCTGATGCTGCTGCTGCGCTG'
	print hmmpred(obs_seq, na3, forw_rerv, 0)
	#obs_seq = 'CTTGTCCGATATGTTCTCTCGCTGTGGCTTCTTGCTCTTGGGGCTAAGTGGATTCTGTATGACTGTCGAACTACATGCAAAAAGT'

	print hmmpred(obs_seq, na3, forw_rerv)


