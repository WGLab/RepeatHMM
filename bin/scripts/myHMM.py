
import string

import numpy as np;
from hmmlearn import hmm

import getTransition_start_emission_prob_2
import getTransition_start_emission_prob_3
import getTransition_start_emission_prob_4
import getTransition_start_emission_prob_5
import getTransition_start_emission_prob_6
import getTransition_start_emission_prob_without0

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


'''
def printHMMmatrix(states, obs_symbols, trainsmat, emisionmat, startprob):
                for si in range(len(states)):
                        if si==0:
                                print ('%6s' % ''),
                                for sj in range(len(states)): print ('%6s' % states[sj]),
                                print ''
                        print ('%6s' % states[si]),
                        sum = 0;
                        for sj in range(len(states)):
                                print ('%.4f' % trainsmat[si][sj]),; sum += trainsmat[si][sj]
                        print ('\t\tsum=%.4f' % sum)
                for si in range(len(states)):
                        if si==0:
                                print ('%6s' % ''),
                                for sj in range(len(obs_symbols)): print ('%6s' % obs_symbols[sj]),
                                print ''
                        print ('%6s' % states[si]),
                        sum = 0;
                        for sj in range(len(obs_symbols)):
                                print ('%.4f' % emisionmat[si][sj]),; sum += emisionmat[si][sj]
                        print ('\t\tsum=%.4f' % sum)
                for si in range(len(states)): print ('%6s' % states[si]),
                print ''
                sum = 0;
                for sj in range(len(states)): print ('%.4f' % startprob[sj]),; sum += startprob[sj]
                print ('\t\tsum=%.4f' % sum)



def getTransition_start_emission_prob_3(na3):
	states = ['N', 'C', 'T', 'G', 'IC', 'IT', 'IG', 'DC', 'DT', 'DG']
	#                  N     C       T       G      IC    IT   IG    DC    DT    DG
	trainsmat = np.array([ 	\
			[0.96, 0.02, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9],          #N 
			[1e-9, 0.0007, 0.8686, 0.0007, 0.11, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #C      0.8566  10.11
			[1e-9, 0.0007, 0.0007, 0.8686, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 0.02],    #T
			[0.02, 0.8486, 0.0007, 0.0007, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9],    #G
			[1e-9, 0.0007, 0.8686, 0.0007, 0.11, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #IC     0.8566   0.11
			[1e-9, 0.0007, 0.0007, 0.8686, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 0.02],    #IT
			[0.02, 0.8486, 0.0007, 0.0007, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9],    #IG
			#[1e-9, 0.0067, 0.9666, 0.0067, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #DC
			#[1e-9, 0.0007, 0.0007, 0.9786, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02],    #DC
			[1e-9, 0.0007, 0.0007, 0.8686, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 0.02],    #DC
			#[1e-9, 0.0067, 0.0067, 0.9666, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02],    #DT
			#[1e-9, 0.9786, 0.0007, 0.0007, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9],    #DT
			[1e-9, 0.8686, 0.0007, 0.0007, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9],    #DT
			#[0.02, 0.9466, 0.0067, 0.0067, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9]]);  #DG
			#[0.02,  0.0007, 0.9586, 0.0007, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9]]);  #DG
			[0.02,  0.0007, 0.8486, 0.0007, 0.11, 1e-9, 1e-9, 1e-9, 0.02, 1e-9]]);  #DG
			#       N    C     T     G     IC    IT   IG     DC    DT    DG 
	startprob = np.array([0.96, 0.02, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9]);

	#                         A      C       G       T	
	emisionmat = np.array([	[0.25, 0.25, 0.25, 0.25],            #N  0
				[0.005, 0.005, 0.005, 0.005],        #C  1
				[0.005, 0.005, 0.005, 0.005],        #T  2
				[0.005, 0.005, 0.005, 0.005],        #G  3
				[0.25, 0.25, 0.25, 0.25],            #IC 4
				[0.25, 0.25, 0.25, 0.25],            #IT 5
				[0.25, 0.25, 0.25, 0.25],            #IG 6
				[0.005, 0.005, 0.005, 0.005],        #DC 7
				[0.005, 0.005, 0.005, 0.005],        #DT 8
				[0.005, 0.005, 0.005, 0.005] ]);     #DG 9
	
	obs_symbols = np.array(['A', 'C', 'G', 'T'])
	for naind in range(len(na3)):
		emind = (np.where(obs_symbols==na3[naind]))[0][0]
		emisionmat[naind+1][emind] = 0.985
		#emisionmat[naind+7][emind] = 1e-9;
		if naind<len(na3)-1:
			afterd = naind + 1;
		else:
			afterd = 0;
		emind = (np.where(obs_symbols==na3[afterd]))[0][0]
		emisionmat[naind+7][emind] = 0.985

	if False:
		for em in emisionmat:
			for ce in em:
				print ('%.4f' % ce),
			print ''

	outputm = False; #True;
	if outputm:
		print 'HMMmatrix1'
		printHMMmatrix(states, obs_symbols, trainsmat, emisionmat, startprob)

	na3 = string.strip(na3);
	state3class = [range(1, len(na3)+1), range(len(na3)+1, 2*len(na3)+1), range(2*len(na3)+1, 3*len(na3)+1)]

	return [trainsmat, startprob, emisionmat, obs_symbols, states, len(states), len(obs_symbols), state3class]

def getTransition_start_emission_prob_3_without0(na3):
	allinfoforhmm = getTransition_start_emission_prob(na3)
	trainsmat = allinfoforhmm[0];
	startprob = allinfoforhmm[1];
	emissionmat = allinfoforhmm[2]
	obs_symbols = allinfoforhmm[3]
	states = allinfoforhmm[4]
	state3class = allinfoforhmm[7]
	
	states = np.delete(states, 0);
	startprob = np.delete(startprob, 0);
	startprob[0] = 0.98

	emissionmat = np.delete(emissionmat, 0, 0);
	
	trainsmat = np.delete(trainsmat, 0, 0);
	trainsmat = np.delete(trainsmat, 0, 1);

	trainsmat[2][0] = 0.8486 + 0.02
	trainsmat[5][0] = 0.8486 + 0.02
	trainsmat[8][1] = 0.8486 + 0.02

	outputm = False; #True;
	if outputm:
		print 'HMMmatrix2'
		printHMMmatrix(states, obs_symbols, trainsmat, emissionmat, startprob)

	return [trainsmat, startprob, emissionmat, obs_symbols, states, len(states), len(obs_symbols), state3class]
'''

def getTransition_start_emission_prob(repPat):
	if len(repPat)==2: return   getTransition_start_emission_prob_2.getTransition_start_emission_prob_2(repPat);
	elif len(repPat)==3: return getTransition_start_emission_prob_3.getTransition_start_emission_prob_3(repPat);
	elif len(repPat)==4: return getTransition_start_emission_prob_4.getTransition_start_emission_prob_4(repPat);
	elif len(repPat)==5: return getTransition_start_emission_prob_5.getTransition_start_emission_prob_5(repPat);
	elif len(repPat)==6: return getTransition_start_emission_prob_6.getTransition_start_emission_prob_6(repPat);
	else: return None 

def getTransition_start_emission_prob_without0(repPat):
	return getTransition_start_emission_prob_without0.getTransition_start_emission_prob_without0(repPat)

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


