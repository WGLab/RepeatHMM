
import os
import sys
import string

import myHMM
import printHMMmatrix
import numpy as np;

len_sensitive_prob = {2:0.8490, 3:0.8486, 4:0.8485, 5:0.8484, 6:0.8480}
outputm = False ; #outputm = True


def getTransition_start_emission_prob_without0(repPat, forprint=False):
        allinfoforhmm = myHMM.getTransition_start_emission_prob(repPat, forprint)
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

        trainsmat[len(repPat)-1][0]   = len_sensitive_prob[len(repPat)] + 0.02
        trainsmat[2*len(repPat)-1][0] = len_sensitive_prob[len(repPat)] + 0.02
        trainsmat[3*len(repPat)-2][0] = len_sensitive_prob[len(repPat)] + 0.02
        trainsmat[3*len(repPat)-1][1] = len_sensitive_prob[len(repPat)] + 0.02

        if outputm or forprint:
                print 'HMMmatrix2'
                printHMMmatrix.printHMMmatrix(states, obs_symbols, trainsmat, emissionmat, startprob)

        return [trainsmat, startprob, emissionmat, obs_symbols, states, len(states), len(obs_symbols), state3class]




