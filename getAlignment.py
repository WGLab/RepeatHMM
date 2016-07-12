
from Bio import pairwise2;
import string

import sys
#sys.path.append('/home/qianliu/project/HTT_CAG_repeat/UnsymmetricPairAlignment')
import UnsymmetricPairAlignment; 
#from UnsymmetricPairAlignment import *

from myheader import *

def getBasePair():
	bp = {}
	bp['A'] = 'T'
	bp['C'] = 'G'
	bp['G'] = 'C'
	bp['T'] = 'A'

	return bp;
def getComplementary(bp, na):
	if not bp.has_key(na.upper()): print na, na.upper(), bp;
	return bp[na.upper()]

def getComplementary3(bp, na3):
	c3 = ''
	for li in range(len(na3)):
		c3 += getComplementary(bp, na3[li]);
	return (c3[::-1]);

def getPattern(na3, forw_rerv):
        if forw_rerv[0]=='-':
                bp = getBasePair()
                return getComplementary3(bp, na3)
	else: return na3;

def getPerfRep(querystr, na3):
	t = int(len(querystr)/2);
	return na3*t


def correctSeq(na3, querystr, forw_rerv, match=2, mismatch=-10, gap=-1, gap_extension=-1):
	na3 = getPattern(na3, forw_rerv)
	perfectstr = getPerfRep(querystr, na3);
	                                                 
	#match:          2;	----->    3	----->	2
	#mismatch:      -1;     ----->   -2	----->	-10
	#gap            -1;     ----->   -2	----->	-1
	#gap extension: -2      ----->   -3	----->	-1
	la = pairwise2.align.localms(perfectstr, querystr, match, mismatch, gap, gap_extension)
	newstr = ''; score=0;
	for a in la:
		p,q,s,b,e=a;
		q = string.strip(q, '-')
		if score<s: score = s; newstr = q;

		#if len(newstr)<len(q): newstr = q;
		#if score<s: score = s; print score, len(newstr);

	return q;

#def myUnsymmetricPairAlignment(na3, querystr, forw_rerv, match=10, mismatch=-9, gap_in_perf=-2, gap_in_read=-13, gap_before_after = -1, bandw=1000, isprint=0):
def myUnsymmetricPairAlignment(na3, querystr, forw_rerv, match=2, mismatch=-2, gap_in_perf=-2, gap_in_read=-12, gap_before_after = -3, bandw=1000, isprint=0):
	na3 = getPattern(na3, forw_rerv)
	perfectstr = getPerfRep(querystr, na3);

	#bp = getBasePair()
	#if forw_rerv[0]=='-':
	#	na3 = getComplementary3(bp, na3)
	#t = int(len(querystr)/2);
	#perfectstr = na3*t;

	newstr = UnsymmetricPairAlignment.correctedByunsymmetricPairWiseAlignment(perfectstr, len(perfectstr), querystr, len(querystr), match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandw,isprint);
	if len(newstr)<len(querystr)/3:	
		print 'In (myUnsymmetricPairAlignment) Warning!!!! too short string from alignment\n', querystr, '\n', newstr
		return querystr
	return newstr;
	
	#newstr = UnsymmetricPairAlignment.unsymmetricPairWiseAlignment(perfectstr, len(perfectstr), querystr, len(querystr), match, mismatch, gap_in_perf, gap_in_read, gap_before_after,isprint);
	
	#if len(newstr.split(';')[2])<len(querystr)/3:
	#	print 'In (myUnsymmetricPairAlignment) Warning!!!! too short string from alignment\n', querystr, '\n', newstr, '\n', newstr.split(';')[2]
	#	return querystr
	#
	#return newstr.split(';')[2]

if __name__=='__main__':
	pat = 'CCG'
	querystr = 'GCCGCCGCGGCGCCCGCCTGCCGCGCCGCTGCCGCGCCCGCAGCAGCCGCTGCCGC'

	pat = 'CTG'
	querystr = 'CTGAGTGTTGAGCTGGCTGACTAGATTGCTTGGTATCGTGCTTTTGCGCTG'
	querystr = 'CATGCTGCTGCTGGCTTCCCGCTGCTGGGTTTTTTTGTTAGTTAATGCTTTTTGCTTGCATGTCTG'
	
	forw_rerv = '+'

	bandw = 50;
	#match=2, mismatch=-10, gap=-1, gap_extension=-1)
	newstr = correctSeq(pat, querystr, '+', match, mismatch, gap_in_perf, gap_in_perf);
	print pat*(int(len(querystr)/2))
	print newstr;
	print querystr


	print ''
	match=10; mismatch=-9; gap_in_perf=-2; gap_in_read=-13; gap_before_after = -1;
	perf = pat*(int(len(querystr)/2))
	print UnsymmetricPairAlignment.unsymmetricPairWiseAlignment(perf, len(perf), querystr, len(querystr), match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandw, 1)


	print '\n', myUnsymmetricPairAlignment(pat, querystr, forw_rerv, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandw, 1)

