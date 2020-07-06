
import os;
import sys;
import string;

#from .myheader import *
from . import myheader

def handletrf_line_sp(line, mdict, otherOpt):
	lsp = line.split()
	#print lsp
	chr, repstr, repend = lsp[:3];
	repstr = int(repstr); repend = int(repend)
	if otherOpt.has_key('chr') and (not otherOpt['chr']==None):
		if not chr==otherOpt['chr']:
			return;
	if otherOpt.has_key('pos') and (not otherOpt['pos']==None):
		if not ((otherOpt['pos'][0]==-1 or repstr>otherOpt['pos'][0]) and (otherOpt['pos'][1]==-1 or repend<otherOpt['pos'][1])):
			return;

	replen, copynum = lsp[4:6];
	repele = lsp[15]

	if not mdict.has_key(chr):
		mdict[chr] = {}
	if not mdict[chr].has_key(repstr):
		mdict[chr][repstr] = {}
	else:
		if mdict[chr][repstr].has_key(repele):
			if myheader.cur_M_STAT <= myheader.M_INFO:
				print ('Warning!!! duplciate repeat', lsp)
				print (chr, repstr, repele, mdict[chr][repstr][repele][:-1])
				print ('                           ', mdict[chr][repstr][repele][-1])

	#									 	 0	    1		   2		  3		  4		    5   6
	mdict[chr][repstr][repele] = [chr, repstr, repend, repele, '+'+copynum, '', lsp]


def readDict(fn, mdict, handle_line, handleoptions=None):
	if not os.path.isfile(fn):
		if myheader.cur_M_STAT <= myheader.M_ERROR: print ('Error !!! no file', fn)
		return ;

	fr = open(fn, 'r')
	line = fr.readline();
	while line:
		if line[0]=='#':
			line = fr.readline(); continue;
 
		line = string.strip(line)

		handle_line(line, mdict, handleoptions)

		line = fr.readline();

def getTRF(moptions):
	if moptions.has_key('bedfile') and (not moptions['bedfile']==None):
		bedfile = moptions['bedfile']
	else:
		if moptions.has_key('hg') and (not moptions['hg']==None):
			hg = moptions['hg']
		else: hg = 'hg38'
		bedfile = moptions['stsBasedFolder']+'/'+hg+'/'+hg+'.trf.bed'
	
	#print 'bedfile', moptions, bedfile
	trfdict = {}
	if os.path.isfile(bedfile):
		readDict(bedfile, trfdict, handletrf_line_sp, moptions)
	else:
		if myheader.cur_M_STAT <= myheader.M_ERROR: 
			print ('No trf.bed file', bedfile)
	#print 'len(trfdict[0])', len(trfdict[trfdict.keys()[0]])
	return trfdict

def handlerepeat_line_sp(line, mdict, otherOpt):
	lsp = line.split(',')
	#  0      1      2          3       4        5       6      7  
	repname, chr, start_pos, end_pos, reppat, strand, range, others = lsp[:8]
	start_pos = int(start_pos)
	end_pos = int(end_pos)

	repname = repname.lower()
	#                  0       1         2        3      4       5      6
	mdict[repname] = [chr, start_pos, end_pos, reppat, strand, range, others, lsp]

def getPredefinedMicrosatellites(moptions):
	if moptions.has_key('pafile') and (not moptions['pafile']==None):
		pafile = moptions['pafile']
	else:
		if moptions.has_key('hg') and (not moptions['hg']==None):
			hg = moptions['hg']
		else: hg = 'hg38'
		pafile = moptions['stsBasedFolder']+'/'+hg+'/'+hg+'.predefined.pa'

	#print 'pafile', pafile, moptions
	trfdict = {}
	if os.path.isfile(pafile):
		readDict(pafile, trfdict, handlerepeat_line_sp, moptions)
	else:
		if myheader.cur_M_STAT <= myheader.M_ERROR: print ('No pa file ', pafile)

	#print 'len(trfdict)', len(trfdict)
	return trfdict


if __name__=='__main__':
	if len(sys.argv)>1:
		hg = sys.argv[1];
	else: hg = 'hg38'
	if len(sys.argv)>2:
		pos_para = sys.argv[2];
	else: pos_para = None;
	
	moptions = {}
	moptions['hg'] = hg
	moptions['stsBasedFolder'] = '../reference_sts'
	moptions['chr'] = 'chr4'
	moptions['pos'] = [3065504, 3086766]
	
	trfdict = getTRF(moptions)
	chrkeys = trfdict.keys(); chrkeys.sort();
	for ck in chrkeys:
		poskeys = trfdict[ck].keys(); poskeys.sort();
		for pk in poskeys:
			elekeys = trfdict[ck][pk].keys(); elekeys.sort();
			for ek in elekeys:
				print (ck, pk, ek, trfdict[ck][pk][ek][:-1])


	print ('')
	mdict = getPredefinedMicrosatellites(moptions)
	rnkeys = mdict.keys(); rnkeys.sort();
	for rk in rnkeys:
		print (rk, mdict[rk][:5]);

	
