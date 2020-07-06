
import os;
import sys;
import string;
import copy
import re;
import logging

#from .myheader import *
from . import myheader
from . import printHMMmatrix
from . import myBAMhandler
from . import myHMM
from . import myGaussianMixtureModel



def replaceSpecificCh(oldc, newc, mstr):
	newlen = len(mstr); oldlen = newlen
	while True:
		mstr = mstr.replace(oldc, newc)
		newlen = len(mstr); 
		if oldlen==newlen: break;
		oldlen = newlen
	return mstr;

def formatOriginalSeqId(oldid):
	newid = replaceSpecificCh('\t', ' ', oldid)
	newid = replaceSpecificCh('  ', ' ', newid)
	newid = newid.replace(' ', '.')

	newid = replaceSpecificCh('__', '_', newid)

	return newid

def getNewFileName(oldfn, uniq_id, mlist, specifiedOptions):
	defaultFolder = specifiedOptions['align']
	if not os.path.isdir(defaultFolder): os.system('mkdir '+defaultFolder)

	curfolder_ind = string.rfind(mlist[0], '/')
	if curfolder_ind==-1: mlist[0] = defaultFolder+mlist[0]
	mlist[0] = defaultFolder+mlist[0][curfolder_ind:]

	if string.find(oldfn, uniq_id)==-1:
		mlist.insert(-1, uniq_id)
	
	cur_final_fname = ''.join(mlist)
	if os.path.isfile(cur_final_fname):
		if myheader.cur_M_STAT <= myheader.M_INFO:  
			print ('Warning!!! filename exist '+cur_final_fname)
			logging.warning('Warning!!! filename exist '+cur_final_fname)
	
	return cur_final_fname;

def obtainFAFromSAM(samfn, uniq_id, specifiedOptions):
	resfn = getNewFileName(samfn, uniq_id, [samfn, '.fa'], specifiedOptions)
	fr = open(samfn, 'r')
	fw = open(resfn, 'w');
	if myheader.cur_M_STAT <= myheader.M_DEBUG: 
		print ('obtainFAFromSAM', samfn, resfn)
		logging.info('obtainFAFromSAM from '+samfn +' to '+ resfn)

	bamseq = {}

	line = fr.readline();
	while line:
		line = string.strip(line);
		lsp = line.split();

		myid = formatOriginalSeqId(lsp[0])
		myseq = lsp[9];

		if bamseq.has_key(myid):
			if len(bamseq[myid])<len(myseq):
				bamseq[myid] = myseq
		else:
			bamseq[myid] = myseq

		line = fr.readline();
	fr.close();

	bamseqkeys = bamseq.keys(); bamseqkeys.sort()
	for bamk in bamseqkeys:
		fw.write('>'+bamk+'\n')
		fw.write(bamseq[bamk]+'\n')
	fw.close();

	if not os.path.isfile(resfn):
		if myheader.cur_M_STAT <= myheader.M_INFO:
			print ('Error!!! There is no result file for '+samfn+'. The result file name is '+resfn)
			logging.info('Error!!! There is no result file for '+samfn+'. The result file name is '+resfn)

	return resfn

def obtainFAFromFQ(fqfn, uniq_id, specifiedOptions):
	resfn = getNewFileName(fqfn, uniq_id, [fqfn, '.fa'], specifiedOptions)
	fr = open(fqfn, 'r')
	fw = open(resfn, 'w');

	if myheader.cur_M_STAT <= myheader.M_DEBUG: 
		print ('obtainFAFromFQ', fqfn, resfn)
		logging.info(''.join(['obtainFAFromFQ from ', fqfn, ' to ', resfn]))
	
	line = fr.readline();
	while line:
		line = string.strip(line);
		if len(line)>0 and line[0]=='@':
			line = ''.join(['>', formatOriginalSeqId(line[1:])]);
			fw.write(line+'\n');
			
			line = string.strip(fr.readline());
			fw.write(line+'\n');
			
			line = fr.readline();
			line = fr.readline();
			
			line = fr.readline();
	
	fr.close();
	fw.close();

	if not os.path.isfile(resfn):
		if myheader.cur_M_STAT <= myheader.M_INFO:
			print ('Error!!! There is no result file for '+fqfn+'. The result file name is '+resfn)
			logging.info(''.join(['Error!!! There is no result file for ', fqfn, '. The result file name is ', resfn]))
	
	return resfn

def getSamePattern1(repPat):
	patset = []
	
	for i in range(len(repPat)):
		if i==0: patset.append(repPat);
		else: patset.append(repPat[i:]+repPat[:i])
	
	mybp = myHMM.getBasePair();
	ptsize = len(patset)
	for i in range(ptsize):
		patset.append(getComplementary(mybp, patset[i]))
	
	return patset

def getSamePattern(moreOptions, commonOptions):
	patset = []
	if commonOptions['CompRep']=='0':
		patset.extend(getSamePattern1(moreOptions['repPat']))
	else: 
		potentialp = []
		for i in range(len(commonOptions['CompRep'])):
			compkeys = commonOptions['CompRep'][i].keys();
			if i==0:
				for ck in compkeys:
					potentialp.append([ck])
			else:
				curlen = len(potentialp)
				for j in range(curlen):
					potentialp[j].append(compkeys[0])
				if len(compkeys)>1:
					for ck in compkeys[1:]:
						for j in range(curlen):
							potentialp.append(copy.deepcopy(potentialp[j]))
							potentialp[-1][-1] = ck
		for curp in potentialp:
			patset.extend(getSamePattern1(''.join(curp)))

	return patset

def containSimilarPattern(commonOptions, patset, curdrep):
	hasFound = False;
	if commonOptions['CompRep']=='0':
		if curdrep in patset: hasFound = True;
	else:
		curPatSize = len(patset[0])
		for exrep in patset:
			if (not curdrep.find(exrep)==-1) and curPatSize/float(len(curdrep))>0.5:
				hasFound = True;
				break;
			if len(curdrep)>=curPatSize:
				identicalNum = 0;
				for i in range(curPatSize):
					if curdrep[i]==exrep[i]: identicalNum += 1
				if identicalNum/float(len(curdrep))>0.5: 
					hasFound = True;
					break;
	return hasFound

def runTRF(fafn, TRFOptions, minRepBWTSize, RepeatTimeThr, moreOptions, commonOptions, specifiedOptions): 
	trf_suf = TRFOptions.replace('_', ' ')
	
	resfn = getNewFileName(fafn, moreOptions['unique_file_id'], [fafn, '_', TRFOptions, '.res'], specifiedOptions)
	#trf_cmd = ''.join(['myTRF ', fafn, ' ', trf_suf, ' -h -ngs > ', resfn])
	trf_cmd = ''.join(['trf ', fafn, ' ', trf_suf, ' -h -ngs > ', resfn])

	patset = getSamePattern(moreOptions, commonOptions)
	if commonOptions['outlog'] <= myheader.M_WARNING:
		if commonOptions['outlog'] <= myheader.M_INFO: print ('patset', patset)
		logging.info('patset ' + str(patset))

	if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2): print (trf_cmd)
	logging.info(trf_cmd)
	os.system(trf_cmd);
	
	repInfo = {}
	
	if not os.path.isfile(resfn):
		if commonOptions['outlog'] <= myheader.M_WARNING:
			print ('Warning!!!! there is no TRF resutls for '+fafn+'. And the result file is '+resfn)
			logging.info(''.join(['Warning!!!! there is no TRF resutls for ',fafn,'. And the result file is ',resfn]))
	else:
		fr = open(resfn,'r')

		curkey = ''
		line = fr.readline();
		while line:
			line = string.strip(line)
			if line[0]=='@':
				curkey = formatOriginalSeqId(string.strip(line[1:]));
				if repInfo.has_key(curkey):
					if commonOptions['outlog'] <= myheader.M_WARNING:
						print ('Warning!!!! Duplicate key: '+curkey)
						logging.info(''.join(['Warning!!!! Duplicate key: ', curkey]))
				else:
					repInfo[curkey] = {}
			else:
				lsp = line.split();
				#[startLoc, endLoc]: 1-index (start from 1) and both tails included
				startLoc = int(lsp[0])
				endLoc = int(lsp[1])
				
				repeatTimes = float(lsp[3]);
				if repeatTimes<RepeatTimeThr: 
					line = fr.readline();
					continue;

				if not containSimilarPattern(commonOptions, patset, lsp[13]):
					line = fr.readline();
					continue;

				curusedkeys = findCloseKeys(repInfo, curkey, minRepBWTSize, startLoc, endLoc)

				if len(curusedkeys)==0:
					repInfo[curkey][startLoc] = [endLoc]
				else:
					curusedkeys.append(startLoc)
					curusedkeys.sort()
					smallestStrk = curusedkeys[0]
					largestEndk = endLoc
					for existk in curusedkeys:
						if repInfo[curkey].has_key(existk) and repInfo[curkey][existk][0]>largestEndk: 
							largestEndk = repInfo[curkey][existk][0]
					
					if not repInfo[curkey].has_key(smallestStrk):
						repInfo[curkey][smallestStrk] = [largestEndk]
					else:
						repInfo[curkey][smallestStrk][0] = largestEndk
					
					for existk in curusedkeys[1:]:
						if repInfo[curkey].has_key(existk):
							del repInfo[curkey][existk]

				if repInfo[curkey].has_key(startLoc):
					if endLoc>repInfo[curkey][startLoc][0]:
						repInfo[curkey][startLoc][0] = endLoc
						
			line = fr.readline();

		fr.close();

		checkRepInfo(repInfo, minRepBWTSize)

		moreOptions['RemList']['trf_resfn'] =(resfn)
		moreOptions['trf_resfn'] =(resfn)
	
	return repInfo

def checkRepInfo(repInfo, minRepBWTSize):
	seqKeys = repInfo.keys();
	for seqk in seqKeys:
		if len(repInfo[seqk])==0: 
			del repInfo[seqk]
			continue;
		posKeys = repInfo[seqk].keys();
		for posk in posKeys:
			existkyes = findCloseKeys(repInfo, seqk, minRepBWTSize, posk, repInfo[seqk][posk][0])
			if len(existkyes)==0 or (posk not in existkyes):
				if myheader.cur_M_STAT <= myheader.M_WARNING:
					print (('Error!!! the current key (%d) is not in the list for %s' % (posk, seqk)), existkyes)
					logging.info(''.join([('Error!!! the current key (%d) is not in the list for %s' % (posk, seqk)), existkyes]))
			if len(existkyes)>1:
				if myheader.cur_M_STAT <= myheader.M_WARNING:
					print ('Waring!!!! two keys are two close for '+seqk+': ', posk, existkyes)
					logging.info(''.join(['Waring!!!! two keys are two close for ',seqk,': ', posk, existkyes]))

def findCloseKeys(repInfo, curkey, minRepBWTSize, startLoc, endLoc):
	startKeys = repInfo[curkey].keys(); startKeys.sort();
	curusedkeys = []
	for stk in startKeys:
		if stk-minRepBWTSize<startLoc<repInfo[curkey][stk][0]+minRepBWTSize or startLoc-minRepBWTSize<stk<endLoc+minRepBWTSize:
			curusedkeys.append(stk)
	return curusedkeys

def splitFA(uniq_id, commonOptions, moreOptions, specifiedOptions):
	fafn = moreOptions['fafile']

	TRFOptions = commonOptions['TRFOptions']
	minRepBWTSize = commonOptions['minRepBWTSize']
	minTailSize = commonOptions['minTailSize']
	RepeatTimeThr = commonOptions['RepeatTime']

	splitfn = getNewFileName(fafn, uniq_id, [fafn, TRFOptions, '.fa'], specifiedOptions)
	spfnbam = getNewFileName(fafn, uniq_id, [fafn, TRFOptions, '_sorted.bam'], specifiedOptions);

	split_norep_fn = getNewFileName(fafn, uniq_id, [fafn, TRFOptions, '_split_norep.fa'], specifiedOptions);
	norep_split_fw = open(split_norep_fn, 'w')
	moreOptions["RemList"]['split_norep_fn'] = split_norep_fn
	moreOptions['split_norep_fn'] = split_norep_fn

	moreOptions['splitfn']=(splitfn)
	moreOptions['spfnbam']=(spfnbam)
	moreOptions['RemList']['splitfn']=(splitfn)
	moreOptions['RemList']['spfnbam']=(spfnbam)

	repInfo = runTRF(fafn, TRFOptions, minRepBWTSize, RepeatTimeThr, moreOptions, commonOptions, specifiedOptions)
	splitfnfw = open(splitfn, 'w')

	splitInfo = {}
	no_repeat_id_list = []
	repeat_id_list = []
	short_repeat_id_list = []
	
	fnfr = open(fafn, 'r');
	seqname = fnfr.readline();  seqcont = ''
	while seqname:
		seqname = formatOriginalSeqId(string.strip(seqname))
		if len(seqname)>0 and seqname[0]=='>':
			haveRep = False;

			seqcont = string.strip(fnfr.readline())
			
			seqk = string.strip(seqname[1:])
			if repInfo.has_key(seqk):
				splitnum = len(repInfo[seqk])+1;
				splitInfo[seqk] = [splitnum, repInfo[seqk], {}, seqcont]

				poskeys = repInfo[seqk].keys(); poskeys.sort();

				start_id = 0; end_id=len(repInfo[seqk])-1
				while start_id<=end_id:
					if start_id==0:
						upstreamsize = poskeys[start_id];
					else:
						upstreamsize = poskeys[start_id] - repInfo[seqk][poskeys[start_id-1]][0];
					if end_id==len(repInfo[seqk])-1:
						downstreamsize = len(seqcont)-repInfo[seqk][poskeys[end_id]][0]
					else:
						downstreamsize = poskeys[end_id+1]-repInfo[seqk][poskeys[end_id]][0]
					if upstreamsize>minTailSize and downstreamsize>minTailSize: break;
					if not upstreamsize>minTailSize: start_id += 1
					if not downstreamsize>minTailSize: end_id -= 1

				repeat_id_list.append([seqk, start_id, end_id, upstreamsize, downstreamsize])

				if upstreamsize>minTailSize and downstreamsize>minTailSize:
					haveRep = True;
					for isp in range(start_id, end_id+1):
						if isp==0:
							pre_end_pos = 0;
						else: 
							pre_end_pos = repInfo[seqk][poskeys[isp-1]][0]
						cur_start_pos = poskeys[isp]

						curseqc = seqcont[pre_end_pos:cur_start_pos]
						cur_id_suf = ''.join([seqname, '__', str(splitnum), '.', str(isp)]);
						if splitInfo[seqk][2].has_key(cur_id_suf):
							if commonOptions['outlog'] <= myheader.M_WARNING and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
								print ('Warning!!! duplciate key in split: '+cur_id_suf)
								logging.info(''.join(['Warning!!! duplciate key in split: ', cur_id_suf]))
						splitInfo[seqk][2][cur_id_suf] = [pre_end_pos, cur_start_pos-1]
						splitfnfw.write(cur_id_suf+'\n')
						splitfnfw.write(curseqc+'\n')
					isp = end_id
					pre_end_pos = repInfo[seqk][poskeys[isp]][0]; 
					if isp<len(repInfo[seqk])-1:
						cur_start_pos = poskeys[isp+1]
						curseqc = seqcont[pre_end_pos:cur_start_pos]
					else:
						cur_start_pos = -1
						curseqc = seqcont[pre_end_pos:]
					cur_id_suf = ''.join([seqname, '__', str(splitnum), '.', str(isp+1)]);
					splitfnfw.write(cur_id_suf+'\n')
					splitfnfw.write(curseqc+'\n')
					if splitInfo[seqk][2].has_key(cur_id_suf):
						if commonOptions['outlog'] <= myheader.M_WARNING:
							print ('Warning!!! duplciate key in split: '+cur_id_suf)
							logging.info(''.join(['Warning!!! duplciate key in splitl: ', cur_id_suf]))
					splitInfo[seqk][2][cur_id_suf] = [pre_end_pos, pre_end_pos+len(curseqc)-1]
				else:
					short_repeat_id_list.append([seqk, upstreamsize, downstreamsize, len(seqcont)])
					if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
						print ('Warning!!! cannot find non-repeat flanking region long enough for ', seqk, repInfo[seqk], upstreamsize, downstreamsize, len(seqcont))
						logging.info(' '.join(['Warning!!! cannot find non-repeat flanking region long enough for ', seqk, str(repInfo[seqk]), str(upstreamsize), str(downstreamsize), str(len(seqcont))]))
					if upstreamsize<0 or downstreamsize<0:
						if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
							print ('Warning!!! negative streamsize: ', upstreamsize, downstreamsize, len(seqcont), repInfo[seqk][poskeys[-1]][0], ' for '+seqname)
							logging.info(' '.join(['Warning!!! negative streamsize: ', str(upstreamsize), str(downstreamsize), str(len(seqcont)), str(repInfo[seqk][poskeys[-1]][0]), ' for '+seqname]))
			else:
				if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
					print ('Info: No repeat information for '+seqname + ' in '+fafn, '>>>', seqcont)
					logging.info(' '.join(['Info: No repeat information for ',seqname , ' in ',fafn, '>>>', seqcont]))

			if not haveRep:
				if seqk in no_repeat_id_list:
					if commonOptions['outlog'] <= myheader.M_WARNING: print ('Warning!!! Duplicate no_repeat_id_list(split)', seqk);
				else:
					no_repeat_id_list.append(seqk)

					norep_split_fw.write(seqname+'\n');
					norep_split_fw.write(seqcont+'\n');

		seqname = fnfr.readline();
	
	moreOptions['no_repeat_id_list'] = no_repeat_id_list
	
	fnfr.close();
	splitfnfw.close();

	if commonOptions['outlog'] <= myheader.M_DEBUG:
		print ('splitfn',splitfn)
		print ('spfnbam',spfnbam)
		#print 'splitInfo',len(splitInfo)
		#splitInfo_keys = splitInfo.keys(); 
		#for sp_k in splitInfo_keys:
		#	print '\t', sp_k, splitInfo[sp_k]
		print ('no_repeat_id_list', len(no_repeat_id_list), '\n', no_repeat_id_list)
		print ('short_repeat_id_list', len(short_repeat_id_list))
		#for short_id in short_repeat_id_list:
		#	print '\t', short_id
		print ('repeat_id_list', len(repeat_id_list))
		#for rep_id in repeat_id_list:
		#	print '\t', rep_id
	
	if (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):	
		logging.info(' '.join(['splitfn',splitfn]))
		logging.info(' '.join(['spfnbam',spfnbam]))
		logging.info(' '.join(['no_repeat_id_list', str(len(no_repeat_id_list))]))
		logging.info(' '.join(['short_repeat_id_list', str(len(short_repeat_id_list))]))
		logging.info(' '.join(['repeat_id_list', str(len(repeat_id_list))]))

	return [splitfn, spfnbam, splitInfo]

def splitFQ(uniq_id, commonOptions, moreOptions, specifiedOptions):
	fqfn = moreOptions['fqfile']
	fafn = obtainFAFromFQ(fqfn, uniq_id, specifiedOptions)
	moreOptions['RemList']['fafile']=fafn
	moreOptions['fafile']=fafn
	res = splitFA(uniq_id, commonOptions, moreOptions, specifiedOptions)

	return res;

def splitSAM(uniq_id, commonOptions, moreOptions, specifiedOptions):
	samfn = moreOptions['samfile']
	fafn = obtainFAFromSAM(samfn, uniq_id, specifiedOptions);
	moreOptions['RemList']['fafile']=fafn
	moreOptions['fafile']=fafn
	
	res = splitFA(uniq_id, commonOptions, moreOptions, specifiedOptions)
	
	return res;

def splitBAM(uniq_id, commonOptions, moreOptions, specifiedOptions):
	chr, startpos, endpos = None, None, None

	if moreOptions.has_key('chr'): chr = moreOptions['chr']
	if moreOptions.has_key('startpos'): startpos = moreOptions['startpos']
	if moreOptions.has_key('endpos'): endpos = moreOptions['endpos']

	bamfn = moreOptions['bamfile']
	bamfn_ind = ''.join([bamfn,'.bai'])
	if not os.path.isfile(bamfn_ind):
		if commonOptions['outlog'] <= myheader.M_WARNING:
			print ('Warning!!!!!!!! the input BAM file not indexed', bamfn)
			logging.info(''.join(['Warning!!!!!!!! the input BAM file not indexed', bamfn]))

	samfn = getNewFileName(bamfn, uniq_id, [bamfn, '.sam'], specifiedOptions)
	if chr==None or startpos==None or endpos==None:
		mview = ''.join(['samtools view ', bamfn, '>', samfn])
	else:
		mview = ''.join(['samtools view ', bamfn, ' ', chr, ':', str(startpos), '-', str(endpos), '>', samfn])

	moreOptions['samfile']=samfn
	moreOptions['RemList']['samfile']=samfn
 
	if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
		print (mview)
	if (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2): logging.info(' '.join(['splitBAM1=', mview]))
	os.system(mview)

	if os.path.getsize(samfn)==0:
		if commonOptions['outlog'] <= myheader.M_WARNING:
			logging.info(mview+'\n')
			logging.info('The file %s have zero size in the first occurrence of splitBAM.\nTry without chr' % samfn)
			#print ('The file %s have zero size\nTry without chr' % samfn)
		if chr==None or startpos==None or endpos==None:
			mview = ''.join(['samtools view ', bamfn, '>', samfn])
		else:
			 mview = ''.join(['samtools view ', bamfn, ' ', chr[3:], ':', str(startpos), '-', str(endpos), '>', samfn])
		if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
			print (mview)
		if (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2): logging.info(' '.join(['splitBAM2=', mview]))
		os.system(mview)

	if os.path.getsize(samfn)==0:
		logging.info('The file %s have zero size: %s in the second occurrence of splitBAM.\n' % (samfn, mview))
		#print ('The file %s have zero size: %s' % (samfn, mview))
		moreOptions[myheader.FATAL_key] = True
		#os.system('rm '+samfn)
		return None

	res = splitSAM(uniq_id, commonOptions, moreOptions, specifiedOptions);
	
	return res;

def reAlign(splitfn, hgfile, splitfn_sorted, moreOptions, commonOptions, specifiedOptions):
	#cmd = (template_bwamem_cmd % (commonOptions['BWAMEMOptions'], hg_reference_and_index, hgfile, splitfn, splitfn_sorted))
	cmd = (myheader.template_bwamem_cmd % (commonOptions['BWAMEMOptions'], '', hgfile, splitfn, splitfn_sorted))
	#print 't1', cmd
	
	if (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
		if commonOptions['outlog'] <= myheader.M_INFO:	print (cmd)
		logging.info(' '.join(['reAlign', cmd]))
	os.system(cmd);

	if not os.path.isfile(splitfn_sorted):
		if commonOptions['outlog'] <= myheader.M_WARNING:
			print ('Warning!!! no file for '+splitfn_sorted)
			logging.info(' '.join(['Warning!!! no file for ', splitfn_sorted]))
	else:
		moreOptions['RemList']['spfnbam'] = splitfn_sorted
		moreOptions['spfnbam'] = splitfn_sorted

	cmd = 'samtools index '+splitfn_sorted
	os.system(cmd)
	if not os.path.isfile(splitfn_sorted+'.bai'):
		if commonOptions['outlog'] <= myheader.M_WARNING:
			print ('Warning!!! no bai file for '+splitfn_sorted)
			logging.info(' '.join(['Warning!!! no bai file for ', splitfn_sorted]))
	else:
		moreOptions['RemList']['spfnbam_bai'] = splitfn_sorted+'.bai'
		moreOptions['spfnbam_bai'] = splitfn_sorted+'.bai'

	if commonOptions['outlog'] <= myheader.M_DEBUG:
		mview = ''.join(['samtools view ', splitfn_sorted, '>', splitfn_sorted, '.sam'])
		print (mview)
		os.system(mview)
		moreOptions['RemList']['spfnbam_sam'] = splitfn_sorted+'.sam'
		moreOptions['spfnbam_sam'] = splitfn_sorted+'.sam'

def getRegioinInBAM(commonOptions, specifiedOptions, moreOptions):
	chr = moreOptions['chr']
	repeatName = moreOptions['repeatName']
	gene_start_end = moreOptions['gene_start_end']
	repeat_start_end = moreOptions['repeat_start_end']
	repPat = string.strip(moreOptions['repPat'])
	forw_rerv = moreOptions['forw_rerv']

	unique_file_id = specifiedOptions['unique_file_id']

	spfnbam = moreOptions['spfnbam']

	alignfolder = specifiedOptions['align']; #'align/'
	if not os.path.isdir(alignfolder): os.system('mkdir '+alignfolder)

	alignfile = alignfolder + repeatName + unique_file_id +'.alignment.sam'
	get_alg_cmd = 'samtools view '+spfnbam+' ' + chr+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+alignfile
	if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2): 
		print (get_alg_cmd)
		logging.info(' '.join(['getRegioinInBAM', get_alg_cmd]))
		logging.info('Running '+get_alg_cmd)
	os.system(get_alg_cmd);
	if os.path.getsize(alignfile)==0:
		if commonOptions['outlog'] <= myheader.M_WARNING and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
			logging.info('The file %s have zero size in the function of getRegioinInBAM.\nTry without chr' % alignfile)
			#print ('The file %s have zero size\nTry without chr' % alignfile)
		get_alg_cmd = 'samtools view '+spfnbam+' ' + chr[3:]+':'+str(gene_start_end[0])+'-'+str(gene_start_end[1])+' > '+alignfile
		if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
			print (get_alg_cmd)
			logging.info('Running '+get_alg_cmd)
		os.system(get_alg_cmd);
	if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
		logging.info('Produced ' + alignfile + ' done!');

	moreOptions['RemList']['spfnbam_sam_interest'] = alignfile
	moreOptions['spfnbam_sam_interest'] = alignfile
	if (not os.path.isfile(alignfile)) or os.path.getsize(alignfile)==0:
		if commonOptions['outlog'] <= myheader.M_INFO:
			logging.error('Cannot produce '+alignfile+' for '+repeatName)
			#print 'Cannot produce '+alignfile+' for '+repeatName
			moreOptions[myheader.FATAL_key] = True

	return alignfile

def getPosOfInterest(alignfile, commonOptions, specifiedOptions, moreOptions):
	alignReader = open(alignfile, 'r')
	chr = moreOptions['chr']

	posDict = {}
	qnameDict = {}

	if not os.path.isfile(alignfile):
		return [posDict, qnameDict]

	cflagDict = {}

	line = alignReader.readline();
	while line:
		line = string.strip(line)

		lsp = line.split('\t')
		qname = lsp[0];
		pos = int(lsp[3])
		cchr = lsp[2];

		if not (chr==cchr or (len(chr)>3 and chr[3:]==cchr) or (len(cchr)>3 and cchr[3:]==chr)):
			continue

		qname_sp = qname.split('__')
		qname_basic = qname_sp[0]
		qname_id = int(qname_sp[1].split('.')[1])

		if not posDict.has_key(qname_basic):
			posDict[qname_basic] = {}
		if not posDict[qname_basic].has_key(qname_id):
			posDict[qname_basic][qname_id] = [{}, qname]
		if (not qnameDict.has_key(qname)) or (len(qnameDict[qname][3])<len(lsp[9])):
			posDict[qname_basic][qname_id][0][(cchr)] = pos

		if not cflagDict.has_key(int(lsp[1])):
			cflagDict[int(lsp[1])] = 0
		cflagDict[int(lsp[1])] += 1

		if not qnameDict.has_key(qname):
			#                    chr          pos  alignInfo     aaInfo
			qnameDict[qname] = [lsp[2], int(lsp[3]), lsp[5],     lsp[9], int(lsp[1])]
		else:
			if myheader.cur_M_STAT <= myheader.M_INFO:
				print ('Duplicate qname: '+qname, line)
				print ('\t\tposDict=', posDict[qname_basic][qname_id])
				print ('\t\tqnameDict=', qnameDict[qname])
				logging.info(' '.join(['Duplicate qname: ', qname, line]))
				logging.info(' '.join(['\t\tposDict=', str(posDict[qname_basic][qname_id])]))
				logging.info(' '.join(['\t\tqnameDict=', str(qnameDict[qname])]))
			if len(qnameDict[qname][3])<len(lsp[9]):
				#                    chr          pos  alignInfo     aaInfo
				qnameDict[qname] = [lsp[2], int(lsp[3]), lsp[5],     lsp[9], int(lsp[1])]

		line = alignReader.readline();
	alignReader.close()

	if myheader.cur_M_STAT <= myheader.M_DEBUG:
		print ('cflagDict')
		cflagDictkeys = cflagDict.keys(); cflagDictkeys.sort()
		for cdk in cflagDictkeys:
			print ('\t', cdk, cflagDict[cdk])

	return [posDict, qnameDict]

def obtain_Flag_group(curf):
	forward_flag = [0, 2048]
	backward_flag = [16, 2064];
	
	if curf in forward_flag: return 1;
	elif curf in backward_flag: return -1
	else: return 0

def getExpectedPosDict(posDict, mchr, qnameDict, splitInfo, moreOptions):
	norep_split_fw = open(moreOptions['split_norep_fn'], 'a')

	expectedPosDict = {}
	qnamekeys = posDict.keys();
	for qk in qnamekeys:
		pos_ordered = {}

		cur_ids = posDict[qk].keys(); cur_ids.sort()
		for cid in cur_ids:
			if posDict[qk][cid][0].has_key(mchr):
					cur_pos = posDict[qk][cid][0][mchr]
					if not pos_ordered.has_key(cur_pos):
						pos_ordered[cur_pos] = [cid]
					else: pos_ordered[cur_pos].append(cid);
		if myheader.cur_M_STAT <= myheader.M_DEBUG: print ('pos_ordered=', pos_ordered, qk, posDict[qk], mchr)
		forward_list = []; forw_list_cflag = []
		backward_list = []; back_list_cflag = []

		pos_keys = pos_ordered.keys(); pos_keys.sort();
		for posk in pos_keys:
			for cid in pos_ordered[posk]:
				find_f = -1;
				curf = obtain_Flag_group(qnameDict[posDict[qk][cid][1]][4])
				for f_ind in range(len(forward_list)):
					fl = forward_list[f_ind]
					if len(fl)>0 and fl[-1][0]+1==cid:
						if myheader.cur_M_STAT <= myheader.M_DEBUG: print ('\t forw', fl, cid, posk, posDict[qk][cid][1], '<', len(qnameDict[posDict[qk][fl[-1][0]][1]][3]), fl[-1][1]+len(qnameDict[posDict[qk][fl[-1][0]][1]][3]))
						if fl[-1][1]<posk:
							existf = forw_list_cflag[f_ind]
							if curf==existf:
								find_f = f_ind
								forward_list[f_ind].append([cid, posk, posDict[qk][cid][1]])
				if find_f==-1:
					forward_list.append([[cid, posk, posDict[qk][cid][1]]])
					forw_list_cflag.append(curf)
				
				find_r = -1;
				for r_ind in range(len(backward_list)):
					rl = backward_list[r_ind]
					if len(rl)>0 and rl[-1][0]-1==cid: 
						if myheader.cur_M_STAT <= myheader.M_DEBUG: print ('\t back', rl, cid, posk, posDict[qk][cid][1], '<', len(qnameDict[posDict[qk][rl[-1][0]][1]][3]), rl[-1][1] + len(qnameDict[posDict[qk][rl[-1][0]][1]][3]))
						if rl[-1][1]<posk:
							existf = back_list_cflag[r_ind]
							if curf==existf:
								find_r = r_ind;
								backward_list[r_ind].append([cid, posk, posDict[qk][cid][1]])
				if find_r==-1:
					backward_list.append([[cid, posk, posDict[qk][cid][1]]])
					back_list_cflag.append(curf)
		
		if myheader.cur_M_STAT <= myheader.M_DEBUG: print ('forward_list=',forward_list)
		if myheader.cur_M_STAT <= myheader.M_DEBUG: print ('backward_list=', backward_list)

		larg_forward_ind=0;
		larg_backward_ind = 0;
		findexp = 0;
		for f_c_ind in range(len(forward_list)):
			if len(forward_list[f_c_ind])>len(forward_list[larg_forward_ind]): 
				larg_forward_ind = f_c_ind
		for r_c_ind in range(len(backward_list)):
			if len(backward_list[r_c_ind])>len(backward_list[larg_backward_ind]):
				larg_backward_ind = r_c_ind
		if len(forward_list)>0 and len(backward_list)>0:
			if len(forward_list[larg_forward_ind])>=len(backward_list[larg_backward_ind]) and len(forward_list[larg_forward_ind])>1:
				expectedPosDict[qk] = [forward_list[larg_forward_ind], 1]; findexp=1
			elif len(forward_list[larg_forward_ind])<len(backward_list[larg_backward_ind]) and len(backward_list[larg_backward_ind])>1:
				expectedPosDict[qk] = [backward_list[larg_backward_ind], -1]; findexp=-1
		elif len(forward_list)>0:
			if len(forward_list[larg_forward_ind])>1: 
				expectedPosDict[qk] = [forward_list[larg_forward_ind], 1]; findexp=1
		elif len(backward_list)>0:
			if len(backward_list[larg_backward_ind])>1: 
				expectedPosDict[qk] = [backward_list[larg_backward_ind], -1]; findexp=-1
		if findexp==0:
			if myheader.cur_M_STAT <= myheader.M_INFO:
				print ('Warning!!! Cannot find expectedPosDict for '+qk, posDict[qk], pos_ordered)
				logging.info(' '.join(['Warning!!! Cannot find expectedPosDict for ', qk, 'posDict=', str(posDict[qk]),'pos_ordered=', str(pos_ordered), 'splitInfo=', str(splitInfo[qk])]))

			if qk in moreOptions['no_repeat_id_list']:
				if myheader.cur_M_STAT <= myheader.M_INFO:
					print ('Warning!!! Duplicate no_repeat_id_list', qk)
					logging.info(' '.join(['Warning!!! Duplicate no_repeat_id_list(pos)', qk]))
			else:
				moreOptions['no_repeat_id_list'].append(qk)
				norep_split_fw.write('>'+qk+'\n');
				norep_split_fw.write(splitInfo[qk][3]+'\n');
		else:
			if myheader.cur_M_STAT <= myheader.M_INFO: print ('expectedPosDict['+qk+']=', expectedPosDict[qk])

	norep_split_fw.close()

	if myheader.cur_M_STAT <= myheader.M_DEBUG:
		print ('expectedPosDict=') 
		expectedPosDict_keys = expectedPosDict.keys();
		for exp_pd_k in expectedPosDict_keys:
			print ('\t', exp_pd_k, expectedPosDict[exp_pd_k])

	return expectedPosDict

def getExpRegionInLongRead(commonOptions, specifiedOptions, moreOptions):
	expRegionInLongRead = {}

	chr = moreOptions['chr']
	repeatName = moreOptions['repeatName']
	gene_start_end = moreOptions['gene_start_end']
	repeat_start_end = moreOptions['repeat_start_end']
	repPat = string.strip(moreOptions['repPat'])
	forw_rerv = moreOptions['forw_rerv']

	repeatFlankLength = commonOptions['repeatFlankLength']

	len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)

	repregion_len_threhold = len_repPat*3; #3;
	if repregion_len_threhold>10: repregion_len_threhold=10
	moreOptions['repeatFlankLength'] = repeatFlankLength

	wrongalign = 0;
	expPosDict = moreOptions['expPosDict']
	qnameDict = moreOptions['qnameDict']
	splitInfo = moreOptions['splitInfo']
	norep_split_fw = open(moreOptions['split_norep_fn'], 'a')

	orig_seq = readRepFASeq(moreOptions)

	expPosDict_keys = expPosDict.keys(); expPosDict_keys.sort()
	minus1num = len(expPosDict_keys);
	for exp_pd_k in expPosDict_keys:
		cur_start_pos_longread = -1;
		cur_end_pos_longread = -1;
		longer = False;
		
		for reg_info in expPosDict[exp_pd_k][0]:
			cur_qname = reg_info[2];
			if not reg_info[1]==qnameDict[cur_qname][1]:
				if commonOptions['outlog'] <= myheader.M_WARNING: 
					print ('Error pos does not match!!! %d=!=%d' % (reg_info[1], qnameDict[cur_qname][1]))
					print ('\t', reg_info, qnameDict[cur_qname][:2], cur_qname)
					logging.info(' '.join([('Error pos does not match!!! %d=!=%d' % (reg_info[1], qnameDict[cur_qname][1]))]))
					logging.info(' '.join(['\t', str(reg_info), str(qnameDict[cur_qname][:2]), cur_qname]))
			
			pos_ref = reg_info[1]
			aligninfo = qnameDict[cur_qname][2]
			aainfo = qnameDict[cur_qname][3]
			
			if (pos_ref > repeat_start_end[0]-repeatFlankLength) and cur_start_pos_longread==-1: 
				wrongalign += 1;
				break;
			
			numreg = re.compile('\d+')
			numinfo = numreg.findall(aligninfo)
			mdireg = re.compile('[MIDNSHPX=]{1}')
			mdiinfo = mdireg.findall(aligninfo)
			
			if not len(numinfo)==len(mdiinfo):
				logging.error('Num is equal to mid' +str(len(numinfo)) + ' '+ str(len(mdiinfo))); continue;
		
			#two tails inclusive	
			cur_split_region_start_end = [splitInfo[exp_pd_k][2]['>'+cur_qname][0], splitInfo[exp_pd_k][2]['>'+cur_qname][1]]
			if expPosDict[exp_pd_k][1]==-1:
				cur_split_region_start_end = [len(orig_seq[exp_pd_k]) - cur_split_region_start_end[1], len(orig_seq[exp_pd_k]) - cur_split_region_start_end[0]]
			cur_split_region_start_end[0] -= 10
			cur_split_region_start_end[1] += 10
			if cur_split_region_start_end[0]<0: 
				cur_split_region_start_end[0] = 0
			if cur_split_region_start_end[1]>len(orig_seq[exp_pd_k]): 
				cur_split_region_start_end[1]=len(orig_seq[exp_pd_k])

			queryind = string.find(orig_seq[exp_pd_k], aainfo, cur_split_region_start_end[0], cur_split_region_start_end[1]);
			if queryind==-1:
				if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
					print ('Warning!!! cannot find tRegion')
					if commonOptions['outlog'] <= myheader.M_DEBUG:
						print ('Warning!!! cannot find tRegion %s in %s = orig_seq[exp_pd_k][%d, %d]' % (aainfo, orig_seq[exp_pd_k][cur_split_region_start_end[0]:cur_split_region_start_end[1]], cur_split_region_start_end[0], cur_split_region_start_end[1]) )
					print (('\texpPosDict[%s]' % exp_pd_k ), expPosDict[exp_pd_k])
					print (('\tsplitInfo[%s][2]' % exp_pd_k), splitInfo[exp_pd_k][2])
					print ('treg_info=', reg_info)
					print (('qnameDict[%s]' % cur_qname), qnameDict[cur_qname][:2])
					print ('\t', string.find(orig_seq[exp_pd_k], aainfo, splitInfo[exp_pd_k][2]['>'+cur_qname][0], splitInfo[exp_pd_k][2]['>'+cur_qname][1]), cur_split_region_start_end, cur_qname, exp_pd_k, string.find(getComplementary(myHMM.getBasePair(), orig_seq[exp_pd_k]), aainfo, splitInfo[exp_pd_k][2]['>'+cur_qname][0], splitInfo[exp_pd_k][2]['>'+cur_qname][1]))
					print ('\taainfo=', aainfo)
					print ('\taainfo=', getComplementary(myHMM.getBasePair(),aainfo))
					print ('\t', orig_seq[exp_pd_k][splitInfo[exp_pd_k][2]['>'+cur_qname][0]:splitInfo[exp_pd_k][2]['>'+cur_qname][1]])
				
					logging.info(' '.join([('Warning!!! cannot find tRegion %s in %s = orig_seq[exp_pd_k][%d, %d]' % (aainfo, orig_seq[exp_pd_k][cur_split_region_start_end[0]:cur_split_region_start_end[1]], cur_split_region_start_end[0], cur_split_region_start_end[1]) )]))
					logging.info(' '.join([('\texpPosDict[%s]' % exp_pd_k ), str(expPosDict[exp_pd_k])]))
					logging.info(' '.join([('\tsplitInfo[%s][2]' % exp_pd_k), str(splitInfo[exp_pd_k][2])]))
					logging.info(' '.join(['reg_info=', str(reg_info)]))
					logging.info(' '.join([('qnameDict[%s]' % cur_qname), str(qnameDict[cur_qname][:2])]))
					logging.info(' '.join(['\t', str(string.find(orig_seq[exp_pd_k], aainfo, splitInfo[exp_pd_k][2]['>'+cur_qname][0], splitInfo[exp_pd_k][2]['>'+cur_qname][1])), str(cur_split_region_start_end), cur_qname, str(exp_pd_k), str(string.find(getComplementary(myHMM.getBasePair(), orig_seq[exp_pd_k]), aainfo, splitInfo[exp_pd_k][2]['>'+cur_qname][0], splitInfo[exp_pd_k][2]['>'+cur_qname][1]))]))
					logging.info(' '.join(['\taainfo=', aainfo]))
					logging.info(' '.join(['\taainfo=', getComplementary(myHMM.getBasePair(),aainfo)]))
					logging.info(' '.join(['\t', str(orig_seq[exp_pd_k][splitInfo[exp_pd_k][2]['>'+cur_qname][0]:splitInfo[exp_pd_k][2]['>'+cur_qname][1]])]))
				continue;

			for n1ind in range(len(numinfo)):
				n1 = int(numinfo[n1ind])
				mdi = mdiinfo[n1ind];
				for n1i in range(n1):
					if mdi=='M' or mdi == '=' or mdi == 'X':
						pos_ref = pos_ref + 1;
						queryind = queryind + 1;
					elif mdi =='I':
						queryind = queryind + 1;
					elif mdi == 'D':
						pos_ref = pos_ref + 1;
					elif mdi == 'S':
						queryind = queryind + 1;
					elif mdi=='H':
						pass; 
					elif mdi=='P':
						pass;
					else:
						if commonOptions['outlog'] <= myheader.M_WARNING:
							logging.warning('Warning unknown CIGAR element ' + str(n1) + ' ' + mdi)
							print ('Warning unknown CIGAR element ' + str(n1) + ' ' + mdi)

					if pos_ref-1 < repeat_start_end[0]-repeatFlankLength:
						cur_start_pos_longread = queryind-1
						if cur_start_pos_longread<0: cur_start_pos_longread = 0
	
					if pos_ref-1 < repeat_start_end[1]+repeatFlankLength:
						cur_end_pos_longread = queryind

		if cur_start_pos_longread==-1 or cur_end_pos_longread==-1 or pos_ref-1<repeat_start_end[1]+repeatFlankLength or pos_ref-1<repeat_start_end[0]-repeatFlankLength:
			longer = False
		if pos_ref-1 >= repeat_start_end[1]+repeatFlankLength:
			longer = True;

		if longer: minus1num -= 1

		addseq = False;
		if longer and cur_end_pos_longread>=0 and cur_start_pos_longread>=0 and cur_end_pos_longread-cur_start_pos_longread >= repregion_len_threhold and cur_end_pos_longread<len(orig_seq[exp_pd_k]):
			currepregion = orig_seq[exp_pd_k][cur_start_pos_longread:cur_end_pos_longread]
			currepregion = currepregion.replace('-', '')
			if len(currepregion)>=repregion_len_threhold:
				addseq = True;
				if (not expRegionInLongRead.has_key(exp_pd_k)):
					expRegionInLongRead[exp_pd_k] = [[longer, currepregion]]
				else:
					expRegionInLongRead[exp_pd_k].append([longer, currepregion])
		if not addseq:	
			if not cur_end_pos_longread<=len(orig_seq[exp_pd_k]):
				if commonOptions['outlog'] <= myheader.M_ERROR:
					print ('Error!!!! cur_end_pos_longread%d greater than %d for %s' % (cur_end_pos_longread, len(orig_seq[exp_pd_k]), exp_pd_k))
			if exp_pd_k in moreOptions['no_repeat_id_list']:
				if commonOptions['outlog'] <= myheader.M_WARNING:
					print ('Warning!!! Duplicate no_repeat_id_list', qk)
					logging.info(' '.join(['Warning!!! Duplicate no_repeat_id_list(longer)', qk]))
			else:
				moreOptions['no_repeat_id_list'].append(exp_pd_k)
				norep_split_fw.write('>'+exp_pd_k+'\n');
				norep_split_fw.write(splitInfo[exp_pd_k][3]+'\n');

	norep_split_fw.close()
	if commonOptions['outlog'] <= myheader.M_INFO:
		print ('Used sp long reads: %d-%d=%d' % (len(expPosDict_keys), minus1num, len(expPosDict_keys)-minus1num))

	if commonOptions['outlog'] <= myheader.M_DEBUG:
		expRegionInLongRead_keys = expRegionInLongRead.keys();
		expRegionInLongRead_keys.sort()
		for epk in expRegionInLongRead_keys:
			print ('ReLongRe=', epk, expRegionInLongRead[epk])

	return [expRegionInLongRead, wrongalign]

def getComplementary(mybp, seq):
	newseq = []
	for li in range(1, len(seq)+1):
		newseq.append(mybp[seq[-li]])
	return ''.join(newseq)

def readRepFASeq(moreOptions):
	expPosDict = moreOptions['expPosDict']

	fn = moreOptions['fafile']

	orig_seq = {}
	fr = open(fn, 'r');
	mybp = myHMM.getBasePair();
	
	line1 = fr.readline();
	while line1:
		line1 = string.strip(line1);
		if len(line1)==0:
			line1 = fr.readline();
			continue;
		curkey = ''
		if (line1[0]=='>'):
			curkey = formatOriginalSeqId(line1[1:])
		if not curkey=='':
			line1 = string.strip(fr.readline());

			if expPosDict.has_key(curkey):
				if orig_seq.has_key(curkey):
					if myheader.cur_M_STAT <= myheader.M_INFO:
						print ('Duplicate seq key '+curkey)
						logging.info(' '.join(['Duplicate seq key ', curkey]))
				else:
					if expPosDict[curkey][1]==-1:
						line1 = getComplementary(mybp, line1)

					orig_seq[curkey] = (line1)

		line1 = fr.readline();

	fr.close;
	return orig_seq

def findRegionOfInterest(commonOptions, specifiedOptions, moreOptions):
	alignfile = getRegioinInBAM(commonOptions, specifiedOptions, moreOptions)
	if moreOptions.has_key(myheader.FATAL_key) and moreOptions[myheader.FATAL_key]==True:
		return None

	posDict, qnameDict = getPosOfInterest(alignfile, commonOptions, specifiedOptions, moreOptions)

	splitInfo = moreOptions['splitInfo']
	norep_split_fw = open(moreOptions['split_norep_fn'], 'a')

	splitInfo_keys = splitInfo.keys();
	for splitk in splitInfo_keys:
		if posDict.has_key(splitk):
			pass;
		else:
			if splitk in moreOptions['no_repeat_id_list']:
				pass
			else:
				if splitk not in moreOptions['no_repeat_id_list']:
					moreOptions['no_repeat_id_list'].append(splitk)
					norep_split_fw.write('>'+splitk+'\n');
					norep_split_fw.write(splitInfo[splitk][3]+'\n');
	norep_split_fw.close()

	if os.path.isfile(alignfile):
		expPosDict = getExpectedPosDict(posDict, moreOptions['chr'], qnameDict, splitInfo, moreOptions)
	
		moreOptions['expPosDict'] = expPosDict
		moreOptions['qnameDict'] = qnameDict
		expRegionInLongRead = getExpRegionInLongRead(commonOptions, specifiedOptions, moreOptions)
		del moreOptions['expPosDict']
		del moreOptions['qnameDict']
	else:
		expRegionInLongRead = [{}, 0]

	return expRegionInLongRead

def findRepeatCountALongRead(commonOptions, specifiedOptions, moreOptions):
	expRegionInLongRead, wrongalign = moreOptions['expRegionInLongRead']
	
	chr = moreOptions['chr']
	repeat_start_end = moreOptions['repeat_start_end']
	repPat = string.strip(moreOptions['repPat'])
	forw_rerv = moreOptions['forw_rerv']
	repeatFlankLength = moreOptions['repeatFlankLength']
	repeatName = moreOptions['repeatName']

	MinSup = commonOptions['MinSup']
	
	len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)

	ref_repeat = (repeat_start_end[1]-repeat_start_end[0]+1)/float(len_repPat)
	hmmoptions = myBAMhandler.getHMMOptions(repeatFlankLength, repPat, forw_rerv, commonOptions)

	rep_predata, rep_mfadata, rep_sufdata = myBAMhandler.getGene(repeatName, chr, repeat_start_end, specifiedOptions['unique_file_id'], specifiedOptions['analysis_file_id'], commonOptions['hgfile'], repeatFlankLength, specifiedOptions)

	commonOptions['rep_flanking_data'] = rep_predata, rep_sufdata

	rptrue = []; rpfalse = []; orignial = [];
	repeatsKeys = expRegionInLongRead.keys(); repeatsKeys.sort()
	for currep_key in repeatsKeys:
		curseqalign = expRegionInLongRead[currep_key]
		for currep in curseqalign: 
			newstr = currep[1]

			pre0 = 0; predstats=''
			if len(newstr)<commonOptions['MaxRep']: #*len_repPat:
				#print 'Realignment', repeat_start_end, chr
				newstr, pre0, predstats = myBAMhandler.getUnsymAlignAndHMM(repPat, forw_rerv, repeatFlankLength, hmmoptions, currep[1], commonOptions, currep_key)
			else: logging.warning('The sequence is too long: '+str(len(newstr))+' '+chr+' '+repeatName+' '+repPat+' '+str(currep[0])+' reads name:'+currep_key+" "+str(commonOptions['MaxRep'])+" "+str(commonOptions['MaxRep']*len_repPat))
			orignial.append([currep[1], pre0, predstats]);
			if currep[0] and commonOptions['outlog'] <= myheader.M_DEBUG: #int(len(newstr)/float(len_repPat)) in [8, 13]: #currep[0] and isTest:
				print (currep_key, ('pred_rep=%.3f' % (len(newstr)/float(len_repPat))))
				if len(newstr)/float(len_repPat)>80: print ('more than 80')
				else: print ('')
				print ('HMM1:') 
				prst = []
				for i in range(pre0): prst.append('-')
				prst.append(newstr)
				print (''.join(prst))
				print ('HMM2:', predstats)
				print ('HMM3:', currep[1])
			
			currep[1] = newstr
			if currep[0]: #
				rptrue.append(len(currep[1])/float(len_repPat)) #
			else:
				rpfalse.append(len(currep[1])/float(len_repPat)) #

	rptrue.sort(); rpfalse.sort()
	trstr = 'true ' + str(len(rptrue)) + ' [';
	for rpt in rptrue:
		trstr = trstr + ('%.0f,' % rpt)
	trstr = trstr[:-1] + ']'
	logging.debug(trstr)

	p2, allocr = myGaussianMixtureModel.get2Peaks(rptrue, MinSup, commonoptions=commonOptions)

	if len(rpfalse)>0:
		flstr = 'fals ' + str(len(rpfalse)) + ' ['
		for rpf in rpfalse:
			flstr = flstr + ('%.0f,' % rpf)
		flstr = flstr[:-1] + ']'
		logging.debug(flstr);

		logging.info('ref_repeat ' + ('%.0f' % ref_repeat) +'\t'+repPat+'\t'+forw_rerv);

	return [repeatName, ref_repeat, p2, allocr, len(rptrue), len(rpfalse)+wrongalign]

def detectRepRegion(commonOptions, specifiedOptions, moreOptions):
	fn = moreOptions['fafqfile']
	fafqtype = moreOptions['fafqtype']

	uniq_id = specifiedOptions['unique_file_id'] 

	TRFOptions = commonOptions['TRFOptions']
	minRepBWTSize = commonOptions['minRepBWTSize']
	minTailSize = commonOptions['minTailSize']
	hgfile = commonOptions['hgfile']
 
	if fafqtype=='fa': 
		moreOptions['fafile'] = fn
		splitfn, spfnbam, splitInfo = splitFA(uniq_id, commonOptions, moreOptions, specifiedOptions) #
	elif fafqtype=='fq':
		moreOptions['fqfile'] = fn
		splitfn, spfnbam, splitInfo = splitFQ(uniq_id, commonOptions, moreOptions, specifiedOptions)
	elif fafqtype=='bam':
		moreOptions['bamfile'] = fn
		mret = splitBAM(uniq_id, commonOptions, moreOptions, specifiedOptions)
		if moreOptions.has_key(myheader.FATAL_key) and moreOptions[myheader.FATAL_key]==True:
			return None
		splitfn, spfnbam, splitInfo = mret
	elif fafqtype=='sam':
		moreOptions['samfile'] = fn
		splitfn, spfnbam, splitInfo = splitSAM(uniq_id, commonOptions, moreOptions, specifiedOptions)
	else: 
		if commonOptions['outlog'] <= myheader.M_ERROR:
			print ('Wrong type of fafaq '+fafqtype)
			logging.info(' '.join(['Wrong type of fafaq ', str(fafqtype)])) 
		return repAllesInfo
		
	reAlign(splitfn, hgfile, spfnbam, moreOptions, commonOptions, specifiedOptions);

	moreOptions['splitInfo'] = splitInfo
	expRegionInLongRead = findRegionOfInterest(commonOptions, specifiedOptions, moreOptions);
	if moreOptions.has_key(myheader.FATAL_key) and moreOptions[myheader.FATAL_key]==True:
		return None
	del moreOptions['splitInfo']

	return expRegionInLongRead

def getNonRepeatAlignment(commonOptions, specifiedOptions, moreOptions):
	if moreOptions['fafqtype'] in ['bam', 'sam']:
		moreOptions['nonRepeat_alignfn_sp'] = moreOptions['samfile']
	else:
		nonRepeat_alignfn_bam = getNewFileName(moreOptions['fafile'], specifiedOptions['unique_file_id'], [moreOptions['fafile'], '_nonRepeat.sam'] , specifiedOptions)
		split_norep_fn = moreOptions['split_norep_fn']
		
		#cmd = (template_bwamem_cmd2 % (hg_reference_and_index, commonOptions['hgfile'], split_norep_fn, nonRepeat_alignfn_bam))
		cmd = (tmyheader.emplate_bwamem_cmd2 % ('', commonOptions['hgfile'], split_norep_fn, nonRepeat_alignfn_bam))
		#print 't2', cmd
		if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
			print (cmd);
			logging.info('getNonRepeatAlignment: '+cmd)
		os.system(cmd)
		moreOptions['RemList']['nonRepeat_alignfn_bam'] = nonRepeat_alignfn_bam
		moreOptions['nonRepeat_alignfn_bam'] = nonRepeat_alignfn_bam

		cmd = 'samtools index '+nonRepeat_alignfn_bam
		os.system(cmd)
		if not os.path.isfile(nonRepeat_alignfn_bam+'.bai'):
			if commonOptions['outlog'] <= myheader.M_WARNING:
				print ('Warning!!! no bai file for '+nonRepeat_alignfn_bam)
				logging.info(' '.join(['Warning!!! no bai file for ', nonRepeat_alignfn_bam]))
		else:
			moreOptions['RemList']['nonRepeat_alignfn_bam_bai'] = nonRepeat_alignfn_bam+'.bai'
			moreOptions['nonRepeat_alignfn_bam_bai'] = nonRepeat_alignfn_bam+'.bai'
	
		gene_start_end = moreOptions['gene_start_end']
		chr = moreOptions['chr']
		nonRepeat_alignfn_bam_sp = ''.join([nonRepeat_alignfn_bam,'_sp.sam'])
		mview = ''.join(['samtools view ', nonRepeat_alignfn_bam, ' ', chr, ':', str(gene_start_end[0]), '-', str(gene_start_end[1]), '>', nonRepeat_alignfn_bam_sp])
		if (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
			if commonOptions['outlog'] <= myheader.M_INFO: print ('getNonRepeatAlignment', mview)
			logging.info(''.join(['getNonRepeatAlignment', mview]))
		os.system(mview)

		moreOptions['RemList']['nonRepeat_alignfn_sp'] = nonRepeat_alignfn_bam_sp
		moreOptions['nonRepeat_alignfn_sp'] = nonRepeat_alignfn_bam_sp

	nonRepeatinLongRead = getNonRepeatinLongRead(commonOptions, specifiedOptions, moreOptions)

	del moreOptions['no_repeat_id_list']
	
	return nonRepeatinLongRead

def getNonRepeatinLongRead(commonOptions, specifiedOptions, moreOptions):
	no_repeat_id_list = moreOptions['no_repeat_id_list']
	nonRepeat_alignfn = moreOptions['nonRepeat_alignfn_sp']

	nonRepeatinLongRead = {}

	chr = moreOptions['chr']
	repeatName = moreOptions['repeatName']
	gene_start_end = moreOptions['gene_start_end']
	repeat_start_end = moreOptions['repeat_start_end']
	repPat = string.strip(moreOptions['repPat'])
	forw_rerv = moreOptions['forw_rerv']

	repeatFlankLength = commonOptions['repeatFlankLength']

	len_repPat = printHMMmatrix.get_len_repPat(repPat, commonOptions)

	repregion_len_threhold = len_repPat*3; #3;
	if repregion_len_threhold>10: repregion_len_threhold=10

	wrongalign = 0;

	non_rep_reader = open(nonRepeat_alignfn, 'r')
	line = non_rep_reader.readline();
	while line:
		line = string.strip(line)

		lsp = line.split('\t')
		exp_pd_k = formatOriginalSeqId(lsp[0])
		
		if exp_pd_k not in  no_repeat_id_list:
			line = non_rep_reader.readline();
			continue
		
		cchr = lsp[2]
		pos_ref = int(lsp[3])
		aligninfo = lsp[5]
		aainfo = lsp[9]
		
		cur_start_pos_longread = -1;
		cur_end_pos_longread = -1;
		longer = False;

		if not (chr==cchr or (len(chr)>3 and chr[3:]==cchr) or (len(cchr)>3 and cchr[3:]==chr)):
			continue

		if not (cchr==chr or cchr==chr[3:]):  
			if commonOptions['outlog'] <= myheader.M_WARNING:
				logging.error('Not same ' + cchr +' ' + chr);
				print ('Error!!! Not same ' + cchr +' ' + chr)
			line = non_rep_reader.readline();
			continue;
	
		if (pos_ref > repeat_start_end[0]-repeatFlankLength) and cur_start_pos_longread==-1:
			if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
				print ('Info: Cannot cover upstream: ', pos_ref, repeat_start_end[0], repeatFlankLength, cur_start_pos_longread)
				logging.info(' '.join(['Info: Cannot cover upstream: ', str(pos_ref), str(repeat_start_end[0]), str(repeatFlankLength), str(cur_start_pos_longread)]))
			wrongalign += 1;
			line = non_rep_reader.readline();
			continue;

		numreg = re.compile('\d+')
		numinfo = numreg.findall(aligninfo)
		mdireg = re.compile('[MIDNSHPX=]{1}')
		mdiinfo = mdireg.findall(aligninfo)
		
		if not len(numinfo)==len(mdiinfo):
			if commonOptions['outlog'] <= myheader.M_INFO:
				logging.error('Num is equal to mid' +str(len(numinfo)) + ' '+ str(len(mdiinfo)));
				print ('Num is equal to mid' +str(len(numinfo)) + ' '+ str(len(mdiinfo)))
			wrongalign += 1;
			line = non_rep_reader.readline();
			continue;

		queryind = 0;
		for n1ind in range(len(numinfo)):
			n1 = int(numinfo[n1ind])
			mdi = mdiinfo[n1ind];
			for n1i in range(n1):
				if mdi=='M' or mdi == '=' or mdi == 'X':
					pos_ref = pos_ref + 1;
					queryind = queryind + 1;
				elif mdi =='I':
					queryind = queryind + 1;
				elif mdi == 'D':
					pos_ref = pos_ref + 1;
				elif mdi == 'S':
					queryind = queryind + 1;
				elif mdi=='H':
					pass;
				elif mdi=='P':
					pass;
				else:
					if commonOptions['outlog'] <= myheader.M_WARNING:
						logging.warning('Warning unknown CIGAR element ' + str(n1) + ' ' + mdi)
						print ('Warning unknown CIGAR element ' + str(n1) + ' ' + mdi)

				if pos_ref-1 < repeat_start_end[0]-repeatFlankLength:
					cur_start_pos_longread = queryind-1;
					if cur_start_pos_longread<0: cur_start_pos_longread = 0
				if pos_ref-1 < repeat_start_end[1]+repeatFlankLength:
					cur_end_pos_longread = queryind

		if cur_start_pos_longread==-1 or cur_end_pos_longread==-1 or pos_ref-1<repeat_start_end[1]+repeatFlankLength or pos_ref-1<repeat_start_end[0]-repeatFlankLength:
			longer = False
		if pos_ref-1 >= repeat_start_end[1]+repeatFlankLength:
			longer = True;
		
		if longer and cur_end_pos_longread>=0 and cur_start_pos_longread>=0 and cur_end_pos_longread-cur_start_pos_longread >= repregion_len_threhold and cur_end_pos_longread<len(aainfo):
			currepregion =  aainfo[cur_start_pos_longread:(cur_end_pos_longread)]
			currepregion = currepregion.replace('-', '')
			if len(currepregion)>= repregion_len_threhold: 
				if not nonRepeatinLongRead.has_key(exp_pd_k):
					nonRepeatinLongRead[exp_pd_k] = [[longer, currepregion]]
				else:
					nonRepeatinLongRead[exp_pd_k].append([longer, currepregion])
			else:
				if commonOptions['outlog'] <= myheader.M_INFO and (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
					print ('Warning!! Non-rep too short ', cur_start_pos_longread, cur_end_pos_longread, len(aainfo), currepregion)
		line = non_rep_reader.readline();
	
	non_rep_reader.close()
	
	return [nonRepeatinLongRead, wrongalign]

def getRepeatCounts(commonOptions, specifiedOptions, moreOptions):
	trf_para_size = len(commonOptions['TRFOptions'].split('_'))
	if trf_para_size<7:
		len_repPat = printHMMmatrix.get_len_repPat(moreOptions['repPat'], commonOptions)
		periodOfinterest = int(len_repPat*2)
		if periodOfinterest - len_repPat>50: periodOfinterest = len_repPat + 50
		commonOptions['TRFOptions'] = commonOptions['TRFOptions']+'_'+str(periodOfinterest)

	if commonOptions['outlog']<=myheader.M_INFO: print ('rep patterns', moreOptions['repPat'], commonOptions['CompRep'])
	logging.info(moreOptions['chr']+' '+str(moreOptions['repeat_start_end']))
	
	moreOptions['startpos'] = moreOptions['gene_start_end'][0]
	moreOptions['endpos'] = moreOptions['gene_start_end'][1]
	if commonOptions['SeqTech']=="Illumina": rep_up_down_size = 300;
	elif commonOptions['SeqTech']=="Pacbio": rep_up_down_size = 5000;
	elif commonOptions['SeqTech']=="Nanopore": rep_up_down_size = 5000;
	if commonOptions['SeqTech'] in ["Illumina", "Pacbio", "Nanopore"]:
		if moreOptions['repeat_start_end'][0]-moreOptions['startpos']>rep_up_down_size: 
			moreOptions['startpos'] = moreOptions['repeat_start_end'][0] - rep_up_down_size
		if moreOptions['endpos']-moreOptions['repeat_start_end'][1]>rep_up_down_size:
			moreOptions['endpos'] = moreOptions['repeat_start_end'][1] + rep_up_down_size
	
	moreOptions['RemList'] = {}
	moreOptions['unique_file_id'] = specifiedOptions['unique_file_id']
	
	expRegionInLongRead = detectRepRegion(commonOptions, specifiedOptions, moreOptions)
	if moreOptions.has_key(myheader.FATAL_key) and moreOptions[myheader.FATAL_key]==True:
		myrm(moreOptions, specifiedOptions, commonOptions);
		return None
	nonRepeatinLongRead = getNonRepeatAlignment(commonOptions, specifiedOptions, moreOptions)

	if commonOptions['outlog'] <= myheader.M_INFO:
		print ('expRegionInLongRead=', str(len(expRegionInLongRead[0])), str(expRegionInLongRead[1]));
		print ('nonRepeatinLongRead=', str(len(nonRepeatinLongRead[0])), str(nonRepeatinLongRead[1]))
	logging.info(' '.join(['expRegionInLongRead=', str(len(expRegionInLongRead[0])), str(expRegionInLongRead[1])]))
	logging.info(' '.join(['nonRepeatinLongRead=', str(len(nonRepeatinLongRead[0])), str(nonRepeatinLongRead[1])]))
	
	expRegionInLongRead[1] += nonRepeatinLongRead[1]
	norepeatkeys = nonRepeatinLongRead[0].keys();
	for nork in norepeatkeys:
		if expRegionInLongRead[0].has_key(nork):
			if commonOptions['outlog'] <= myheader.M_INFO:
				print ('Duplicate key in non-repeat:', nork)
				logging.info(''.join(['Duplicate key:', nork]));
		else:
			expRegionInLongRead[0][nork] = nonRepeatinLongRead[0][nork]

	moreOptions['expRegionInLongRead'] = expRegionInLongRead
	repcount = findRepeatCountALongRead(commonOptions, specifiedOptions, moreOptions)
	del moreOptions['expRegionInLongRead']
	
	myrm(moreOptions, specifiedOptions, commonOptions);

	return repcount

def myrm(moreOptions, specifiedOptions, commonOptions):
	remkeys = moreOptions['RemList'].keys(); 
	for remk in remkeys:
		rem_cmd = ''.join(['rm ', moreOptions['RemList'][remk]])
		if os.path.isfile(moreOptions['RemList'][remk]):
			os.system(rem_cmd)
			if (not specifiedOptions.has_key('thread') or specifiedOptions['thread']<2):
				if commonOptions['outlog'] <= myheader.M_INFO: print (remk, rem_cmd)
				logging.info(' '.join(['rm', remk, rem_cmd]))
		del moreOptions[remk]
	del moreOptions['RemList']
	

if __name__=='__main__':
	moreOptions = {}
	moreOptions['fafqfile'] = '/home/qianliu/project/HTT_CAG_repeat/myTestAlign/test.fq'
	moreOptions['fafqfile'] = '/home/qianliu/project/HTT_CAG_repeat/myTestAlign/sam015.raw.fastq'
	moreOptions['fafqfile'] = '/home/qianliu/project/HTT_CAG_repeat/myTestAlign/sam015.ccs.fastq'
	moreOptions['fafqtype'] = 'fq'

	specifiedOptions = {}

	isfree = False;
	#isfree = True;
	if isfree:
		moreOptions['fafqfile'] = '/home/qianliu/project/HTT_CAG_repeat/myTestAlign/freeze4-all-merge.sort.bam'
		moreOptions['fafqtype'] = 'bam'
		specifiedOptions['unique_file_id'] = '.freeze4-all-merge.sort.bam.test.'
		specifiedOptions['analysis_file_id'] = '.freeze4-all-merge.sort.bam.test.'
	else:
		specifiedOptions['unique_file_id'] = '.test.'
		specifiedOptions['analysis_file_id'] = '.test.'
		
	#moreOptions['chr'] = 'chr14'
	#moreOptions['startpos'] = 92070888
	#moreOptions['endpos'] = 92072403

	commonOptions = {}
	commonOptions['TRFOptions'] = '2_7_7_80_10_50_500'
	commonOptions['TRFOptions'] = '2_7_4_80_10_100_500'
	commonOptions['minRepBWTSize'] = 70 # 10
	commonOptions['minTailSize'] = 70 #30 #70 #10
	commonOptions['hgfile'] = 'hg38.fa'

	moreOptions['chr'] = 'chr14'
	moreOptions['repeatName'] = 'atxn3'

	moreOptions['gene_start_end'] = [92038652, 92106610]

	moreOptions['gene_start_end'] = [92070888, 92072403]
	moreOptions['repeat_start_end'] = [92071011, 92071052]
	moreOptions['repPat'] = 'CAG'
	moreOptions['forw_rerv'] = '-14'

	commonOptions['isGapCorrection'] = 1
	commonOptions['repeatFlankLength'] = 30
	commonOptions['MinSup'] = 3
	commonOptions['RepeatTime'] = 5


	ishtt = False;
	#ishtt = True;
	if ishtt:
		moreOptions['chr'] = 'chr4'
		moreOptions['repeatName'] = 'htt'
		moreOptions['gene_start_end'] = [3074510, 3243960]
		moreOptions['repeat_start_end'] = [3074877, 3074939]
		moreOptions['repPat'] = 'CAG'
		moreOptions['forw_rerv'] = '+19'
		#moreOptions['startpos'] = 3074510
		#moreOptions['endpos'] = 3243960


	repcount = getRepeatCounts(commonOptions, specifiedOptions, moreOptions)
	#detectRepRegion(commonOptions, specifiedOptions, moreOptions)
	
	print (repcount)


