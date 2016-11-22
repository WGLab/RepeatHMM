
import os;
import sys;

import string

import numpy as np;

import getTransition_start_emission_prob_without0
import printHMMmatrix

def getTransition_start_emission_prob_3(repPat):
	avgsub = 0.0007

	repPat = string.strip(repPat);
	if not len(repPat)==3: print 'Error, not 3-acid repeat', repPat;

	states = ['N', 'C', 'T', 'G', 'IC', 'IT', 'IG', 'DC', 'DT', 'DG']
	#                  N     C       T       G      IC    IT   IG    DC    DT    DG
	trainsmat = np.array([ 	\
			[0.96, 0.02,   1e-9,   1e-9,   1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9],          #N 
			[1e-9, avgsub, 0.8686, avgsub, 0.11, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #C
			[1e-9, avgsub, avgsub, 0.8686, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 0.02],    #T
			[0.02, 0.8486, avgsub, avgsub, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9],    #G
			[1e-9, avgsub, 0.8686, avgsub, 0.11, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #IC
			[1e-9, avgsub, avgsub, 0.8686, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 0.02],    #IT
			[0.02, 0.8486, avgsub, avgsub, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9],    #IG
			[1e-9, avgsub, avgsub, 0.8686, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 0.02],    #DC
			[0.02, 0.8486, avgsub, avgsub, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9],    #DT
			[0.02, avgsub, 0.8486, avgsub, 0.11, 1e-9, 1e-9, 1e-9, 0.02, 1e-9]]);  #DG
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
	for naind in range(len(repPat)):
		emind = (np.where(obs_symbols==repPat[naind]))[0][0]
		emisionmat[naind+1][emind] = 0.985
		#emisionmat[naind+7][emind] = 1e-9;
		if naind<len(repPat)-1:
			afterd = naind + 1;
		else:
			afterd = 0;
		emind = (np.where(obs_symbols==repPat[afterd]))[0][0]
		emisionmat[naind+1+len(repPat)*2][emind] = 0.985

	if False:
		for em in emisionmat:
			for ce in em:
				print ('%.4f' % ce),
			print ''

	if getTransition_start_emission_prob_without0.outputm:
		print 'HMMmatrix1'
		printHMMmatrix.printHMMmatrix(states, obs_symbols, trainsmat, emisionmat, startprob)

	repPat = string.strip(repPat);
	state3class = [range(1, len(repPat)+1), range(len(repPat)+1, 2*len(repPat)+1), range(2*len(repPat)+1, 3*len(repPat)+1)]

	return [trainsmat, startprob, emisionmat, obs_symbols, states, len(states), len(obs_symbols), state3class]


if __name__=='__main__':
	getTransition_start_emission_prob_without0.getTransition_start_emission_prob_without0('CAG')

