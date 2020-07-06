
import os;
import sys;

import string
import logging
import numpy as np;

from . import printHMMmatrix

#from .myheader import *
from . import myheader

def getMajor(commonOptions):
	repplist = []
	for i in range(len(commonOptions['CompRep'])):
		compkeys = commonOptions['CompRep'][i].keys(); compkeys.sort()
		repplist.append(compkeys[0])
		for ck in compkeys:
			if commonOptions['CompRep'][i][ck]>commonOptions['CompRep'][i][repplist[-1]]:
				repplist[-1] = ck
	return ''.join(repplist)


def produce_tolerate_mismatch(repPat, commonOptions):
	mismatchstr=[]; bfaf = 2;
	if not commonOptions['CompRep']=='0': 
		if repPat=='':
			repPat = getMajor(commonOptions)

		for i in range(len(commonOptions['CompRep'])):
			compkeys = commonOptions['CompRep'][i].keys(); compkeys.sort()
			if len(compkeys)>1:
				for ck in compkeys:
					if ck==repPat[i]: continue;
					
					start_ind = 0; end_ind = 0;
					bf_1 = True; af_1 = True;
					for j in range(1, bfaf+1):
						if i-j>=0 and bf_1:
							start_ind = -j;
							#if len(commonOptions['CompRep'][i-j])==1:
							#	start_ind = -j;
							#else: bf_1 = False;
						if i+j<len(repPat) and af_1:
							end_ind = j;
							#if len(commonOptions['CompRep'][i+j])==1:
							#	end_ind = j;
							#else: af_1 = False;
					if start_ind==0 and end_ind==0: 
						if commonOptions['outlog'] <= myheader.M_WARNING: print ('Warning!!! both 0 {} {} {} {}'.format( repPat, commonOptions, i, ck))
					else:
						rpstr = ''
						tlstrlist = ['']
						for j in range(start_ind, end_ind+1):
							rpstr = rpstr + repPat[i+j]
							if len(commonOptions['CompRep'][i+j])==1:
								for k in range(len(tlstrlist)):
									tlstrlist[k] = tlstrlist[k]+commonOptions['CompRep'][i+j].keys()[0];
							elif len(commonOptions['CompRep'][i+j])>1:
								newtlstrlist = []
								addcomkeys = commonOptions['CompRep'][i+j].keys();
								for ack in addcomkeys:
									for k in range(len(tlstrlist)):
										if j==0 and (not ack==ck): continue;
										newtlstrlist.append(tlstrlist[k]+ack)
								tlstrlist = newtlstrlist
							else:
								if commonOptions['outlog'] <= myheader.M_WARNING: print ('Fatal error zero comprep {} {} {} {}'.format( j, i+j, commonOptions['CompRep'][i+j], commonOptions['CompRep']))
						for atl in tlstrlist:
							mismatchstr.append(('%s:%s:%d:%d;' % (rpstr, atl, start_ind, end_ind)))

	if commonOptions.has_key('Tolerate_mismatch') and (not commonOptions['Tolerate_mismatch']==None):
		addstr = commonOptions['Tolerate_mismatch'].split(';')
		for ads in addstr:
			mismatchstr.append(ads+';')

	mismnum = len(mismatchstr)
	retmis = ''.join(mismatchstr)[:-1]
	
	return [retmis, len(retmis), mismnum]

def getTransition_start_emission_prob_x(repPat, commonOptions, forprint=False):
	repPat = string.strip(repPat);
	if (commonOptions['CompRep']=='0' and len(repPat)<1) and (len(commonOptions['CompRep'])<1 and (not commonOptions['CompRep']=='0')): return None

	len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)
	if commonOptions['CompRep']=='0': CompRepPat = printHMMmatrix.getCompRepFromSimple(repPat)
	else: CompRepPat = commonOptions['CompRep']

	tol_info = produce_tolerate_mismatch(repPat, commonOptions)
	if commonOptions['outlog'] <= myheader.M_INFO: print ('tol_info={}'.format( tol_info))
	logging.info('tol_info=' + str(tol_info))

	avgsub = 0.0005
	avgsub = myheader.hmm_random_rep_transit/len_repPat
	avgsub = 1e-9
	
	typeOfRepEle = ['', 'I', 'D'];
	repEle = [];
	for rp_ind in range(len_repPat):
		repEle.append(''.join(['r', str(rp_ind+1)]));
	states = ['N'];
	for typRE in typeOfRepEle:
		for rp in repEle:
			states.append(''.join([typRE, rp]));

	if commonOptions.has_key('transitionm') and (not commonOptions['transitionm']==None):
		trainsmat = commonOptions['transitionm']
	else:
		trainsmat = np.full((len(states), len(states)), 1e-9);
		#for N to N
		trainsmat[0][0] = 0.96;
		#for N to rep;
		if not len_repPat<2:
			trainsmat[0][1] = 0.02;
		else: trainsmat[0][1] = 0.04;
		if not len_repPat<2:
			trainsmat[0][1+len(repEle)*2] = 0.02;
		#for rep to N;
		trainsmat[len(repEle)][0] = 0.02;
		trainsmat[len(repEle)*2][0] = 0.02;
		if not len_repPat<2:
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
				trainsmat[len(repEle)*typ_ind+j+1][jind+add_index] = commonOptions['hmm_insert_rate'] #0.11
		#for deletion
		add_index = len(repEle)*2+1;
		for typ_ind in range(len(typeOfRepEle)):
			for j in range(len(repEle)):
				for k in range(1, len(repEle)):
					if typ_ind<len(typeOfRepEle)-1:
						jind = j+k
					else:
						if k>=len(repEle)-1: continue
						jind = j+k+1
					if jind > len(repEle)-1: jind -= len(repEle);
					trainsmat[len(repEle)*typ_ind+j+1][jind+add_index] = commonOptions['hmm_del_rate']**k # 
					if trainsmat[len(repEle)*typ_ind+j+1][jind+add_index]<1e-9: 
						trainsmat[len(repEle)*typ_ind+j+1][jind+add_index] = 1e-9
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
	if len_repPat<2:
		startprob[1] = 0.04
	else:
		startprob[1] = 0.02
		startprob[1+len(repEle)*2] = 0.02
	startprob = np.array(startprob)

	if commonOptions.has_key('emissionm') and (not commonOptions['emissionm']==None):
		emisionmat = commonOptions['emissionm']
	else:
		#emisionmat = np.full((len(repEle)*len(typeOfRepEle)+1, 4), commonOptions['hmm_sub_rate']/4)
		emisionmat = np.full((len(repEle)*len(typeOfRepEle)+1, 5), commonOptions['hmm_sub_rate']/4)
		randrow = [0]
		for j in range(len(repEle)):
			randrow.append(j+len(repEle)+1);
		if len_repPat<2: randrow.append(len(repEle)*len(typeOfRepEle))
		#print randrow
		for rdr in randrow:
			for jcol in range(4):
				#emisionmat[rdr][jcol] = 0.25;
				if not rdr==0:
					emisionmat[rdr][jcol] = 0.25;
				else: emisionmat[rdr][jcol] = 0.2
		for nset in range(len(repEle)*len(typeOfRepEle)+1):
			if not nset==0: emisionmat[nset][4] = 1e-9
			else: emisionmat[nset][4] = 0.2
	
		obs_symbols = np.array(['A', 'C', 'G', 'T', 'N'])
		for naind in range(len_repPat):
			CompRepPatkeys1 = CompRepPat[naind].keys();
			for k1 in CompRepPatkeys1:
				emind = (np.where(obs_symbols==k1))[0][0]
				emisionmat[naind+1][emind] += (1-commonOptions['hmm_sub_rate'])*CompRepPat[naind][k1]
			if len_repPat<2: continue;
			if naind<len_repPat-1:
				afterd = naind + 1;
			else:
				afterd = 0;
			CompRepPatkeys2 = CompRepPat[afterd].keys();
			for k2 in CompRepPatkeys2:
				emind = (np.where(obs_symbols==k2))[0][0]
				emisionmat[naind+1+len_repPat*2][emind] += (1-commonOptions['hmm_sub_rate'])*CompRepPat[afterd][k2]

	if forprint:
		if commonOptions['outlog'] <= myheader.M_INFO: 
			print ('HMMmatrix1')
			printHMMmatrix.printHMMmatrix(states, obs_symbols, trainsmat, emisionmat, startprob)

	state3class = [range(1, len_repPat+1), range(len_repPat+1, 2*len_repPat+1), range(2*len_repPat+1, 3*len_repPat+1)]

	#           0         1           2           3            4                  5             6               7          8   
	return [trainsmat, startprob, emisionmat, obs_symbols, np.array(states), len(states), len(obs_symbols), state3class, tol_info]

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
			print ('\t <%d> 1=%s, 2=%s' % (i, a1[i], a2[i]))
	return isSame
def CompareTwoNumpyArrays(matnames, matindex, fun, matorg, matx):
	for mid in range(len(matnames)):
		a1 = matorg[matindex[mid]];  a2 = matx[matindex[mid]]
		print matnames[mid],
		if fun(a1, a2):
			print (': Same')

def compareMat(matorg, matx):
	matnames = ['trainsmat', 'emisionmat']
	matindex = [    0,          2]
	CompareTwoNumpyArrays(matnames, matindex, compareTwoNumpyArray, matorg, matx)

	matnames = ['startprob', 'obs_symbols']
	matindex = [   1,            3,       ]
	CompareTwoNumpyArrays(matnames, matindex, compareTwoNumpyArray1, matorg, matx)

	CompareTwoNumpyArrays(['states'], [4], compareTwoMat1, matorg, matx)

	if matorg[5]==matx[5]: print ('state_num is same')
	else: print ('state_num is not same: %d, %d' % (matorg[5], matx[5]))
	if matorg[6]==matx[6]: print ('state_num is same')
	else: print ('symbols_num is not same: %d, %d' % (matorg[6], matx[6]))
	
	print ('state3class')
	if compareTwoMat(matorg[7], matx[7]):
		print ('\t: Same')
	

if __name__=='__main__':
	wouldprint = False; 

	commonOptions = {}
	commonOptions['CompRep'] = '0'
	commonOptions['hmm_sub_rate'] = 0.02
	commonOptions['hmm_del_rate'] = 0.02
	commonOptions['hmm_insert_rate'] = 0.11
	
	print ('len(Pattern) = 4')
	matx4   = getTransition_start_emission_prob_x('TATC', commonOptions, wouldprint);
	print ('\nlen(Pattern) = 3')
	matx3   = getTransition_start_emission_prob_x('CAG', commonOptions, wouldprint);
	print ('\nlen(Pattern) = 2')
	matx2   = getTransition_start_emission_prob_x('CG', commonOptions, wouldprint);
	
	print ('\nlen(Pattern) = 5')
	matx5   = getTransition_start_emission_prob_x('TATCG', commonOptions, wouldprint);

	print ('\nlen(Pattern) = 6')
	matx6   = getTransition_start_emission_prob_x('TATCGG', commonOptions, wouldprint);

	wouldprint = True; wouldprint = False;
	
	matx1   = getTransition_start_emission_prob_x('G', commonOptions, wouldprint);

	commonOptions['CompRep'] = 'AlTlT50/C50lClT/C'
	commonOptions['CompRep'] = printHMMmatrix.getCompRep(commonOptions['CompRep'])
	matxy = getTransition_start_emission_prob_x('ATTCT', commonOptions, wouldprint)
	
	print ('\n\ncomplex')
	matxy = getTransition_start_emission_prob_x('', commonOptions, wouldprint)

	commonOptions['CompRep'] = '0'
	matx1   = getTransition_start_emission_prob_x('CG', commonOptions, wouldprint);
	matx6   = getTransition_start_emission_prob_x('TATCGG', commonOptions, wouldprint)
	print (matx6[0][-1])


	wouldprint = True; wouldprint = False; 
	commonOptions['CompRep'] = 'AlTlT40/C50/A10lClT/C'
	commonOptions['CompRep'] = printHMMmatrix.getCompRep(commonOptions['CompRep'])
	matxy = getTransition_start_emission_prob_x('', commonOptions, wouldprint)
	print (matxy[-1])

	commonOptions['CompRep'] = 'ClAlA40/G60'
	commonOptions['CompRep'] = printHMMmatrix.getCompRep(commonOptions['CompRep'])
	matxy = getTransition_start_emission_prob_x('', commonOptions, wouldprint)
	print (matxy[-1])

