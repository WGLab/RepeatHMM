
import os;
import sys;

import string

import numpy as np;

import getTransition_start_emission_prob_without0
import printHMMmatrix

def getTransition_start_emission_prob_6(repPat, forprint=False):
        avgsub = 0.0004

        repPat = string.strip(repPat);
        if not len(repPat)==6: print 'Error, not 6-acid repeat', repPat;

        states = ['N', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'Ir1', 'Ir2', 'Ir3', 'Ir4', 'Ir5', 'Ir6', 'Dr1', 'Dr2', 'Dr3', 'Dr4', 'Dr5', 'Dr5']
        #  N    r1     r2       r3      r4     r5      r6      Ir1   Ir2   Ir3   Ir4   Ir5   Ir6   Dr1   Dr2   Dr3   Dr4   Dr5  Dr6
        trainsmat = np.array([  \
        [0.96, 0.02,   1e-9,   1e-9,   1e-9,   1e-9,   1e-9,   1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9],    #N
        [1e-9, avgsub, 0.8680, avgsub, avgsub, avgsub, avgsub, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9, 1e-9],    #r1
        [1e-9, avgsub, avgsub, 0.8680, avgsub, avgsub, avgsub, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9],    #r2
        [1e-9, avgsub, avgsub, avgsub, 0.8680, avgsub, avgsub, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9],    #r3
        [1e-9, avgsub, avgsub, avgsub, avgsub, 0.8680, avgsub, 1e-9, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #r4
        [1e-9, avgsub, avgsub, avgsub, avgsub, avgsub, 0.8680, 1e-9, 1e-9, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02],    #r5
        [0.02, 0.8480, avgsub, avgsub, avgsub, avgsub, avgsub, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9],    #r6
        [1e-9, avgsub, 0.8680, avgsub, avgsub, avgsub, avgsub, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9, 1e-9],    #Ir1
        [1e-9, avgsub, avgsub, 0.8680, avgsub, avgsub, avgsub, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9],    #Ir2
        [1e-9, avgsub, avgsub, avgsub, 0.8680, avgsub, avgsub, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9],    #Ir3
        [1e-9, avgsub, avgsub, avgsub, avgsub, 0.8680, avgsub, 1e-9, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #Ir4
        [1e-9, avgsub, avgsub, avgsub, avgsub, avgsub, 0.8680, 1e-9, 1e-9, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02],    #Ir5
        [0.02, 0.8480, avgsub, avgsub, avgsub, avgsub, avgsub, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9],    #Ir6
        [1e-9, avgsub, avgsub, 0.8680, avgsub, avgsub, avgsub, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9],    #Dr1
        [1e-9, avgsub, avgsub, avgsub, 0.8680, avgsub, avgsub, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9],    #Dr2
        [1e-9, avgsub, avgsub, avgsub, avgsub, 0.8680, avgsub, 1e-9, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #Dr3
        [1e-9, avgsub, avgsub, avgsub, avgsub, avgsub, 0.8680, 1e-9, 1e-9, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02],    #Dr4
        [0.02, 0.8480, avgsub, avgsub, avgsub, avgsub, avgsub, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9],    #Dr5
        [0.02, avgsub, 0.8480, avgsub, avgsub, avgsub, avgsub, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9, 1e-9]]);  #Dr6
         #                     N     r1    r2    r3    r4   r5    r6    Ir1   Ir2   Ir3   Ir4   Ir5   Ir6   Dr1   Dr2   Dr3   Dr4   Dr5   Dr6
        startprob = np.array([0.96, 0.02, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9]);

        #                         A      C       G       T
        emisionmat = np.array([ [0.25,  0.25,  0.25,  0.25],         #N  0
                                [0.005, 0.005, 0.005, 0.005],        #r1  1
                                [0.005, 0.005, 0.005, 0.005],        #r2  2
                                [0.005, 0.005, 0.005, 0.005],        #r3  3
                                [0.005, 0.005, 0.005, 0.005],        #r4  4
                                [0.005, 0.005, 0.005, 0.005],        #r5  5
                                [0.005, 0.005, 0.005, 0.005],        #r6  6
                                [0.25,  0.25,  0.25,  0.25],         #Ir1 7
                                [0.25,  0.25,  0.25,  0.25],         #Ir2 8
                                [0.25,  0.25,  0.25,  0.25],         #Ir3 9
                                [0.25,  0.25,  0.25,  0.25],         #Ir4 10
                                [0.25,  0.25,  0.25,  0.25],         #Ir5 11
                                [0.25,  0.25,  0.25,  0.25],         #Ir6 12
                                [0.005, 0.005, 0.005, 0.005],        #Dr1 13
                                [0.005, 0.005, 0.005, 0.005],        #Dr2 14
                                [0.005, 0.005, 0.005, 0.005],        #Dr3 15
                                [0.005, 0.005, 0.005, 0.005],        #Dr4 16
                                [0.005, 0.005, 0.005, 0.005],        #Dr5 17
                                [0.005, 0.005, 0.005, 0.005] ]);     #Dr6 18

        obs_symbols = np.array(['A', 'C', 'G', 'T'])
        for naind in range(len(repPat)):
                emind = (np.where(obs_symbols==repPat[naind]))[0][0]
                emisionmat[naind+1][emind] = 0.985
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

        return [trainsmat, startprob, emisionmat, obs_symbols, states, len(states), len(obs_symbols), state3class]



if __name__=='__main__':
	getTransition_start_emission_prob_without0.getTransition_start_emission_prob_without0('TATCGG')

