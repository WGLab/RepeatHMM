
import os;
import sys;

import string

import numpy as np;

import getTransition_start_emission_prob_without0
import printHMMmatrix

def getTransition_start_emission_prob_4(repPat):
        avgsub = 0.0005

        repPat = string.strip(repPat);
        if not len(repPat)==4: print 'Error, not 4-acid repeat', repPat;

        states = ['N', 'r1', 'r2', 'r3', 'r4', 'Ir1', 'Ir2', 'Ir3', 'Ir4', 'Dr1', 'Dr2', 'Dr3', 'Dr4']
        #                  N    r1     r2       r3      r4     Ir1   Ir2   Ir3   Ir4   Dr1   Dr2   Dr3   Dr4
        trainsmat = np.array([  \
                        [0.96, 0.02,   1e-9,   1e-9,   1e-9,   1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9],    #N
                        [1e-9, avgsub, 0.8685, avgsub, avgsub, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9],    #r1
                        [1e-9, avgsub, avgsub, 0.8685, avgsub, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #r2
                        [1e-9, avgsub, avgsub, avgsub, 0.8685, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 0.02],    #r3
                        [0.02, 0.8485, avgsub, avgsub, avgsub, 1e-9, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9, 1e-9],    #r4
                        [1e-9, avgsub, 0.8685, avgsub, avgsub, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9],    #Ir1
                        [1e-9, avgsub, avgsub, 0.8685, avgsub, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #Ir2
                        [1e-9, avgsub, avgsub, avgsub, 0.8685, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 0.02],    #Ir3
                        [0.02, 0.8485, avgsub, avgsub, avgsub, 1e-9, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9, 1e-9],    #Ir4
                        [1e-9, avgsub, avgsub, 0.8685, avgsub, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9],    #Dr1
                        [1e-9, avgsub, avgsub, avgsub, 0.8685, 1e-9, 1e-9, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 0.02],    #Dr2
                        [0.02, 0.8485, avgsub, avgsub, avgsub, 1e-9, 1e-9, 1e-9, 0.11, 0.02, 1e-9, 1e-9, 1e-9],    #Dr3
                        [0.02, avgsub, 0.8485, avgsub, avgsub, 0.11, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9]]);  #Dr4
         #                     N     r1    r2    r3    r4    Ir1  Ir2   Ir3   Ir4   Dr1   Dr2   Dr3   Dr4
        startprob = np.array([0.96, 0.02, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 0.02, 1e-9, 1e-9, 1e-9]);

        #                         A      C       G       T
        emisionmat = np.array([ [0.25,  0.25,  0.25,  0.25],         #N  0
                                [0.005, 0.005, 0.005, 0.005],        #r1  1
                                [0.005, 0.005, 0.005, 0.005],        #r2  2
                                [0.005, 0.005, 0.005, 0.005],        #r3  3
                                [0.005, 0.005, 0.005, 0.005],        #r4  4
                                [0.25,  0.25,  0.25,  0.25],         #Ir1 5
                                [0.25,  0.25,  0.25,  0.25],         #Ir2 6
                                [0.25,  0.25,  0.25,  0.25],         #Ir3 7
                                [0.25,  0.25,  0.25,  0.25],         #Ir4 8
                                [0.005, 0.005, 0.005, 0.005],        #Dr1 9
                                [0.005, 0.005, 0.005, 0.005],        #Dr2 10
                                [0.005, 0.005, 0.005, 0.005],        #Dr3 11
                                [0.005, 0.005, 0.005, 0.005] ]);     #Dr4 12

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

        if getTransition_start_emission_prob_without0.outputm:
                print 'HMMmatrix1'
                printHMMmatrix.printHMMmatrix(states, obs_symbols, trainsmat, emisionmat, startprob)

        state3class = [range(1, len(repPat)+1), range(len(repPat)+1, 2*len(repPat)+1), range(2*len(repPat)+1, 3*len(repPat)+1)]

        return [trainsmat, startprob, emisionmat, obs_symbols, states, len(states), len(obs_symbols), state3class]




if __name__=='__main__':
	getTransition_start_emission_prob_without0.getTransition_start_emission_prob_without0('TATC')

