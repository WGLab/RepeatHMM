
import os;
import sys;
import string;

import time
import resource
import logging


from myheader import *

import myPredefinedPatternReader
import myBAMhandler
import myHMM
import myRepeatReAlignment

import myCommonFun

import multiprocessing


def getAllMicrosatellites(commonOptions, specifiedOptions):
	moptions = {}
	moptions['hg']  = commonOptions['hg']
	moptions['stsBasedFolder'] = commonOptions['stsBasedFolder']
	moptions['outlog'] = commonOptions['outlog']

	scan_region = specifiedOptions['scan_region']
	if not scan_region==None:
		scan_region_sp = scan_region.split(':')
		if len(scan_region_sp)>0:
			moptions['chr'] = scan_region_sp[0]
		else: moptions['chr'] = None
		if len(scan_region_sp)>1:
			pos = scan_region_sp[1].split('-')
			if string.strip(pos[0])=='': pos[0] = -1
			else:	pos[0] = int(pos[0]); 
			if string.strip(pos[1])=='': pos[1] = -1
			else: pos[1] = int(pos[1])
			moptions['pos'] = pos
		else: moptions['pos'] = None
	else:
		moptions['chr'] = None
		moptions['pos'] = None

	#print 'scan_region', scan_region, moptions, commonOptions['Patternfile']
	specifiedOptions['microsatellites'] = [[], []]
	if not commonOptions['Patternfile'] == None:
		for fi in commonOptions['Patternfile']:
			if fi[-3:]=='.pa':
				moptions['pafile'] = fi
				specifiedOptions['microsatellites'][0].append(myPredefinedPatternReader.getPredefinedMicrosatellites(moptions))
			elif fi[-4:]=='.bed':
				moptions['bedfile'] = fi
				specifiedOptions['microsatellites'][1].append(myPredefinedPatternReader.getTRF(moptions))
	else:
		specifiedOptions['microsatellites'][0].append(myPredefinedPatternReader.getPredefinedMicrosatellites(moptions))
		specifiedOptions['microsatellites'][1].append(myPredefinedPatternReader.getTRF(moptions))

	return moptions

def filterMicrosatellites(commonOptions, specifiedOptions, moreOptions):
	allmicro = specifiedOptions['microsatellites']
	moptions = moreOptions['moptions']
	if commonOptions['outlog'] == M_DEBUG: print 'moptions in filterMicrosatellites', moptions
	max_unit_len = 10;

	micronum = 0; #mtotal = 0
	for ti in range(len(allmicro)):
		for mi in range(len(allmicro[ti])):
			if ti==0:
				repkeys = allmicro[ti][mi].keys(); repkeys.sort()
				for repk in repkeys:
					if moptions.has_key('chr') and (not moptions['chr']==None) and (not moptions['chr']==allmicro[ti][mi][repk][0]):
						#print '1.chr', repk, allmicro[ti][mi][repk], moptions['chr'], micronum
						del allmicro[ti][mi][repk]
					elif moptions.has_key('pos') and (not moptions['pos']==None) and (\
						(moptions['pos'][0]>-1 and allmicro[ti][mi][repk][1]<moptions['pos'][0]) or \
						(moptions['pos'][1]>-1 and allmicro[ti][mi][repk][2]>moptions['pos'][1])):
						#print '1.pos', repk, allmicro[ti][mi][repk], moptions['pos'], micronum
						del allmicro[ti][mi][repk]
					elif len(allmicro[ti][mi][repk][3])>max_unit_len:
						#print '1.ele', repk, allmicro[ti][mi][repk], max_unit_len, micronum
						del allmicro[ti][mi][repk]
					else: 
						#print 'keep', repk, allmicro[ti][mi][repk], max_unit_len, moptions['pos'], micronum
						micronum += 1
			else:
				chrkeys = allmicro[ti][mi].keys(); chrkeys.sort()
				for chrk in chrkeys:
					#mtotal = mtotal + len(allmicro[ti][mi][chrk])
					#print ti, mi, chrk, moptions.has_key('chr'), (not moptions['chr']==None), (not moptions['chr']==chrk), moptions['chr']
					if moptions.has_key('chr') and (not moptions['chr']==None) and (not moptions['chr']==chrk):
						if commonOptions['outlog']<=M_WARNING:
							print 'Warning chr del', chrk, chrkeys
						del allmicro[ti][mi][chrk]
						continue;
					startkeys = allmicro[ti][mi][chrk].keys(); startkeys.sort()
					for sk in startkeys:
						if moptions.has_key('pos') and (not moptions['pos']==None) and \
							(moptions['pos'][0]>-1 and sk<moptions['pos'][0]):
							if commonOptions['outlog']<=M_WARNING:
								print 'Warning start-pos del', chrk, sk, allmicro[ti][mi][chrk][sk]
							del allmicro[ti][mi][chrk][sk];
							continue;
						elekeys = allmicro[ti][mi][chrk][sk].keys();
						for ek in elekeys:
							if len(ek)>max_unit_len: 
								del allmicro[ti][mi][chrk][sk][ek];
								continue
							if moptions.has_key('pos') and (not moptions['pos']==None) and \
								(moptions['pos'][1]>-1 and allmicro[ti][mi][chrk][sk][ek][2]>moptions['pos'][1]):

								if commonOptions['outlog']<=M_WARNING:
									print 'Warning end-pos del', chrk, sk, ek, allmicro[ti][mi][chrk][sk][ek]
								del allmicro[ti][mi][chrk][sk][ek]
							else: micronum += 1
						if len(allmicro[ti][mi][chrk][sk])==0:
							del allmicro[ti][mi][chrk][sk]
					if len(allmicro[ti][mi][chrk])==0: del allmicro[ti][mi][chrk]
	#print 'mtotal', mtotal
	return micronum

def printProgress(handlei, micronum, start_time):
	if handlei%50==0 or handlei>micronum-1:
		perc = handlei/float(micronum)*100
		now_time = time.time();
		used_time = now_time-start_time
		tot_time = used_time/perc*100
		
		print ('running time=%.0f (%d%%=%d/%d). Total estimation: %d' % (used_time, int(perc), handlei, micronum, tot_time)); 
		sys.stdout.flush()	

def remove_finished(specifiedOptions):
	maxind = {};
	if not specifiedOptions['continue']==0:
		allmicro = specifiedOptions['microsatellites']
		allfs = os.listdir(''.join([specifiedOptions['scanresfolder'], 'res_', specifiedOptions['analysis_file_id'], '*.txt','.done']));
		mres = {}; mdetail = {}

		for curf in allfs:
			curid = curf[(cuf.rindex(specifiedOptions['analysis_file_id'])+len(specifiedOptions['analysis_file_id'])):-9];
			curidsp = curid.split('_');  #'predef_'
			if not maxind.has_key(curidsp[0]): maxind[curidsp[0]] = -1
			finishedid = int(curidsp[1])
			if finishedid>maxind[curidsp[0]]: maxind[curidsp[0]] = finishedid

			mres, mdetail = myReadScanResults(specifiedOptions, mres, mdetail, curid)
		
		for ti in range(len(allmicro)):
			for mi in range(len(allmicro[ti])):
				if ti==0:
					repkeys = allmicro[ti][mi].keys(); repkeys.sort()
					for rk in repkeys:
						if mres.has_key(rk): del allmicro[ti][mi][rk]
				else:
					chrkeys = allmicro[ti][mi].keys(); chrkeys.sort()
					for chrk in chrkeys:
						startkeys = allmicro[ti][mi][chrk].keys(); startkeys.sort()
						for start_pos in startkeys:
							elekeys = allmicroallmicro[ti][mi][chrk][start_pos].keys();
							haskey = False;
							for repele in elekeys:
								end_pos = str(allmicroallmicro[ti][mi][chrk][start_pos][repele][2]);
								currepkey = ''.join([chrk, ':', start_pos, ':', end_pos, ':', repele]);
								if mres.has_key(currepkey):
									haskey = True;
									break;
							if haskey: del allmicroallmicro[ti][mi][chrk][start_pos]
					if len(allmicroallmicro[ti][mi][chrk])==0: del allmicroallmicro[ti][mi][chrk]
	return maxind
def scan_multiprocess(commonOptions, specifiedOptions):
	multiprocessing.log_to_stderr(logging.DEBUG)

	moptions = getAllMicrosatellites(commonOptions, specifiedOptions)
	moreOptions = {}
	moreOptions['moptions'] = moptions
	micronum = filterMicrosatellites(commonOptions, specifiedOptions, moreOptions)
	allmicro = specifiedOptions['microsatellites']
	maxind = remove_finished(specifiedOptions)
	print 'maxind', maxind
	del specifiedOptions['microsatellites']
	
	avergnum = micronum/specifiedOptions['thread']
	if avergnum>1000: avergnum = 1000;
	minfostr = ('Average microsatellites per thread: %d=%d/%d' % (avergnum, micronum, specifiedOptions['thread']))
	print minfostr
	logging.info(minfostr)

	partitiondict = {};
		
	for ti in range(len(allmicro)):
		for mi in range(len(allmicro[ti])):
			if ti==0:
				repkeys = allmicro[ti][mi].keys(); repkeys.sort()
				splitTask(repkeys, avergnum, partitiondict, 'predef_', maxind, ti, mi, )
			else:
				chrkeys = allmicro[ti][mi].keys(); chrkeys.sort()
				for chrk in chrkeys:
					startkeys = allmicro[ti][mi][chrk].keys(); startkeys.sort()
					splitTask(startkeys, avergnum, partitiondict, chrk+'_', maxind, ti, mi, chrk)

	print 'split taks done', len(partitiondict); sys.stdout.flush()
	unfinishedjobs, finishedjobs, finishedjob_times = distribute_jobs(commonOptions, specifiedOptions, moreOptions, partitiondict, allmicro)
	print 'unfinishedjobs', len(unfinishedjobs), unfinishedjobs; sys.stdout.flush()
	if len(unfinishedjobs)>0:
		moreunfinishedjobs, morefinishedjobs, finishedjob_times = handle_unfinished_jobs(commonOptions, specifiedOptions, moreOptions, unfinishedjobs, partitiondict, avergnum, allmicro, finishedjob_times)
		finishedjobs.update(morefinishedjobs)
		moreunfinishedjobskeys = moreunfinishedjobs.keys()
		if len(moreunfinishedjobskeys)>0:
			print 'Unfinished:!!!!!!!!!!!!!!!!!!!'
			for unf in moreunfinishedjobskeys:
				print '\t', unf, moreunfinishedjobs[unf]
	
	mres = {}; mdetail = {}
	fjkeys = finishedjobs.keys(); fjkeys.sort()
	for fjk in fjkeys:
		mres, mdetail = myCommonFun.myReadScanResults(specifiedOptions, mres, mdetail, fjk);

	return [mres,mdetail]
	

def distribute_jobs(commonOptions, specifiedOptions, moreOptions, partitiondict, allmicro, finishedjob_times={}):
	partitiondictkeys = partitiondict.keys(); partitiondictkeys.sort();
	jobs = {}; usetimes = {}
	finishedjobs = {};
	unfinishedjobs = {}
	if not finishedjob_times.has_key('N90'): finishedjob_times['N90'] = -1;
	if not finishedjob_times.has_key('alltimes'): finishedjob_times['alltimes'] = []
	if not finishedjob_times.has_key('CalNum'): finishedjob_times['CalNum'] = 0
	start_time = time.time();
	pk_ind = 0;

	caltimethre = 25; minrepeatinjobs = 100; msp = 10; pre_usedtime = 0;
	preprint = 0;
	#while pk_ind<len(partitiondictkeys):
	#	while True:
	while True:
			if checkJobs(jobs, usetimes, finishedjobs, finishedjob_times, unfinishedjobs, partitiondict, specifiedOptions) and (not pk_ind<len(partitiondictkeys)): break;

			if len(jobs)<specifiedOptions['thread']:
				#print 'distribute_jobs', len(finishedjob_times['alltimes']), caltimethre, finishedjob_times['N90']
				if len(finishedjob_times['alltimes'])>caltimethre and finishedjob_times['CalNum']<len(finishedjob_times['alltimes']):
					finishedjob_times['alltimes'].sort();
					finishedjob_times['N90'] = finishedjob_times['alltimes'][int(len(finishedjob_times['alltimes'])*0.9)]
					finishedjob_times['CalNum'] = len(finishedjob_times['alltimes']);
				while pk_ind<len(partitiondictkeys) and len(jobs)<specifiedOptions['thread']:
					pk = partitiondictkeys[pk_ind]
					print 'dis', ''.join([specifiedOptions['scanresfolder'], 'res_', specifiedOptions['analysis_file_id'], pk, '.txt','.done'])
					if os.path.isfile(''.join([specifiedOptions['scanresfolder'], 'res_', specifiedOptions['analysis_file_id'], pk, '.txt','.done'])) and specifiedOptions['continue']==0:
						os.system(''.join(['rm ', specifiedOptions['scanresfolder'], 'res_', specifiedOptions['analysis_file_id'], pk, '.txt','.done']))
					moreOptions['curPartition'] = partitiondict[pk]
					ti, mi, chrk = moreOptions['curPartition'][-1]
					if ti==0: specifiedOptions['curmicrosatellites'] = allmicro[ti][mi]
					else: specifiedOptions['curmicrosatellites'] = allmicro[ti][mi][chrk]
					p = multiprocessing.Process(target=scan_part, args=(commonOptions, specifiedOptions, moreOptions, pk,))
					jobs[pk] = p
					usetimes[pk] = [time.time(),time.time()]
					p.start()
					pk_ind = pk_ind + 1

				if len(partitiondict[partitiondictkeys[0]])>minrepeatinjobs:
					cur_time = time.time();
					used_time = cur_time - start_time;
					if len(finishedjobs)>0 and (len(finishedjobs)/msp>preprint or used_time-pre_usedtime>1000):
						pre_usedtime = used_time
						perc = len(finishedjobs)/float(len(partitiondictkeys))
						tot_time = used_time/perc
						memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
						print ('running time=%.0f (%d%%=%d/%d). Total estimation: %d. Mem=%d' % (used_time, int(perc*100), len(finishedjobs), len(partitiondictkeys), tot_time, memres));
						preprint = len(finishedjobs)/msp
						sys.stdout.flush()

			print 'jobs', len(jobs), specifiedOptions['thread'], 'times:', len(finishedjob_times['alltimes']), finishedjob_times['N90']; sys.stdout.flush()
			time.sleep(600);

	#while True:
	#	if checkJobs(jobs, usetimes, finishedjobs, finishedjob_times, unfinishedjobs, partitiondict, specifiedOptions): break;
	#	else: time.sleep(60);
	
	return [unfinishedjobs, finishedjobs, finishedjob_times]

def handle_unfinished_jobs(commonOptions, specifiedOptions, moreOptions, old_unfinishedjobskeys, old_partitiondict, old_avergnum, allmicro, finishedjob_times):
	avergnum = old_avergnum/10;
	if avergnum<1: avergnum = 1
	partiiondict = {}
	for cur_unf_key in old_unfinishedjobskeys:
		if len(old_partitiondict[cur_unf_key])>2:
			ti, mi, chr = old_partitiondict[cur_unf_key][-1]
			splitTask(old_partitiondict[cur_unf_key][:-1], avergnum, partitiondict, cur_unf_key+'_', ti, mi, chr)
	unfinishedjobs, finishedjobs, finishedjob_times = distribute_jobs(commonOptions, specifiedOptions, moreOptions, partitiondict, allmicro, finishedjob_times)
	if len(unfinishedjobs)>0:
		moreunfinishedjobs, morefinishedjobs, finishedjob_times = handle_unfinished_jobs(commonOptions, specifiedOptions, moreOptions, unfinishedjobskeys, partitiondict, avergnum, allmicro, finishedjob_times)
		finishedjobs.update(morefinishedjobs)
		unfinishedjobs = moreunfinishedjobs
	return [unfinishedjobs, finishedjobs, finishedjob_times]

def checkJobs(jobs, usetimes, finishedjobs, finishedjob_times, unfinishedjobs, partitiondict, specifiedOptions):
	jobskeys = jobs.keys(); jobskeys.sort();

	for jk in jobskeys:
		usetimes[jk][1] = time.time()
		if os.path.isfile(''.join([specifiedOptions['scanresfolder'], 'res_', specifiedOptions['analysis_file_id'], jk, '.txt','.done'])):
			print 'Finished: ', ''.join([specifiedOptions['scanresfolder'], 'res_', specifiedOptions['analysis_file_id'], jk, '.txt','.done']), len(jobs), finishedjob_times['N90']; 
			finishedjobs[jk] = 1
			finishedjob_times['alltimes'].append( usetimes[jk][1]- usetimes[jk][0])
			p = jobs[jk]; 
			if p.is_alive(): time.sleep(20);
			p.terminate(); del p; del jobs[jk]; 
			del usetimes[jk];
		else:
			cur_used_time = usetimes[jk][1]- usetimes[jk][0]
			#print jobs[jk].name, cur_used_time, finishedjob_times['N90'], len(finishedjob_times['alltimes']), ';', 
			if finishedjob_times['N90']>0 and cur_used_time>finishedjob_times['N90']*1.5:
				print 'Try to kill', cur_used_time, finishedjob_times['N90']*1.5; sys.stdout.flush()
				jobs[jk].terminate();
				while jobs[jk].is_alive():
					time.sleep(20);
				if not jobs[jk].is_alive():
					unfinishedjobs[jk] = partitiondict[jk]
					p = jobs[jk]; p.is_alive(); p.terminate(); del p; del jobs[jk];
					del usetimes[jk];
				print 'kill'; sys.stdout.flush()
	#print '\tcheckJobs', len(jobs)
	sys.stdout.flush()

	jobskeys = jobs.keys();
	if len(jobskeys)==0 and len(unfinishedjobs)==0: return True;
	else: return False;

def scan_part(commonOptions, specifiedOptions, moreOptions, pk):
	old_analysis_file_id = specifiedOptions['analysis_file_id']
	oldunique_file_id = specifiedOptions['unique_file_id']
	specifiedOptions['analysis_file_id'] = specifiedOptions['analysis_file_id'] + pk
	specifiedOptions['unique_file_id'] = specifiedOptions['unique_file_id'] + pk

	curPartition = moreOptions['curPartition']
	ti, mi, chrk = curPartition[-1]
	allmicro = specifiedOptions['curmicrosatellites']
	mres = {}; mdetail = {}
	moreOptions['mres'] = mres
	moreOptions['mdetail'] = mdetail

	if ti==0:
		for repk in curPartition[:-1]:	
			moreOptions['repeatName'] = repk.lower();
			moreOptions['mgloc'] = allmicro[repk]
			moreOptions['ids'] = False;
			detectRepCounts(commonOptions, specifiedOptions, moreOptions);
	else:
		for sk in curPartition[:-1]:
			elekeys = allmicro[sk].keys();
			for ek in elekeys:
				moreOptions['repeatName'] = ek
				moreOptions['mgloc'] = allmicro[sk][ek]
				moreOptions['ids'] = True;
				detectRepCounts(commonOptions, specifiedOptions, moreOptions);
	specifiedOptions['analysis_file_id'] = old_analysis_file_id
	specifiedOptions['unique_file_id'] = oldunique_file_id
	if specifiedOptions['continue']==0:
		myCommonFun.myWriteScanResults(specifiedOptions, moreOptions['mres'], moreOptions['mdetail'], pk)
	else: myCommonFun.myWriteScanResults(specifiedOptions, moreOptions['mres'], moreOptions['mdetail'], pk, 'a')

def splitTask(repkeys, avergnum, partitiondict, part_name_prefix, maxind, ti, mi, chr=None):
	if avergnum<1: avergnum = 1;
	old_pos = 0; curpart = 0;
	if len(part_name_prefix)>0 and maxind.has_key(part_name_prefix[:-1]):
		curpart = maxind[part_name_prefix[:-1]]+1
	while old_pos<len(repkeys):
		cur_part_name = part_name_prefix+('%09d' % curpart)
		if old_pos+avergnum*1.5<len(repkeys):
			partitiondict[cur_part_name] = repkeys[old_pos:(old_pos+avergnum)]
		else:
			partitiondict[cur_part_name] = repkeys[old_pos:]
		partitiondict[cur_part_name].append([ti, mi, chr])
		curpart = curpart + 1
		old_pos = old_pos + avergnum


def scan(commonOptions, specifiedOptions):
	moptions = getAllMicrosatellites(commonOptions, specifiedOptions)

	mres = {}; mdetail = {}
	moreOptions = {}
	moreOptions['mres'] = mres
	moreOptions['mdetail'] = mdetail
	moreOptions['moptions'] = moptions

	micronum = filterMicrosatellites(commonOptions, specifiedOptions, moreOptions)
	allmicro = specifiedOptions['microsatellites']

	print 'Total size', micronum, allmicro[1][0].keys()

	handlei = 0; start_time = time.time();
	for ti in range(len(allmicro)):
		for mi in range(len(allmicro[ti])):
			if ti==0:
				repkeys = allmicro[ti][mi].keys(); repkeys.sort()
				for repk in repkeys:
					handlei += 1
					moreOptions['repeatName'] = repk.lower();
					moreOptions['mgloc'] = allmicro[ti][mi][repk]
					moreOptions['ids'] = False;
					detectRepCounts(commonOptions, specifiedOptions, moreOptions);
					printProgress(handlei, micronum, start_time)
			else:
				chrkeys = allmicro[ti][mi].keys(); chrkeys.sort()

				for chrk in chrkeys:
					startkeys = allmicro[ti][mi][chrk].keys(); startkeys.sort()
					for sk in startkeys:
						elekeys = allmicro[ti][mi][chrk][sk].keys();
						for ek in elekeys:
							handlei += 1
							moreOptions['repeatName'] = ek
							moreOptions['mgloc'] = allmicro[ti][mi][chrk][sk][ek]
							moreOptions['ids'] = True;
							detectRepCounts(commonOptions, specifiedOptions, moreOptions);
							printProgress(handlei, micronum, start_time)

	return [moreOptions['mres'], moreOptions['mdetail']]

def addSumForAGene(p2, mstr, ind, commonOptions, specifiedOptions, moreOptions):
	myBAMhandler.fixsize2(p2, ind)

	repinfo1 = moreOptions['mgloc']
	chr = repinfo1[0];
	start_pos = str(repinfo1[1]);
	end_pos = str(repinfo1[2]);
	repele = str(repinfo1[3]);

	if moreOptions['ids']:
		currepkey = ''.join([chr, ':', start_pos, ':', end_pos, ':', repele]);
	else:
		currepkey = moreOptions['repeatName']
	detail = ''.join(['', ('%10s' % mstr), ' ', str(p2), '><']);
	retstr = ''.join([str(p2[ind][0]), '/', str(p2[ind][1])]);

	if not moreOptions['mres'].has_key(currepkey):
		moreOptions['mres'][currepkey] = []
	moreOptions['mres'][currepkey].append(retstr)

	if not moreOptions['mdetail'].has_key(currepkey):
		moreOptions['mdetail'][currepkey] = []
	moreOptions['mdetail'][currepkey].append(detail)
	
	#print currepkey, retstr, detail

def detectRepCounts(commonOptions, specifiedOptions, moreOptions):
	retoptions = myBAMhandler.get_Loc1(moreOptions['mgloc'], commonOptions)
	if commonOptions['outlog'] <= M_INFO: print 'mgloc', moreOptions['mgloc']
	#print moreOptions['mgloc'][:5]

	moreOptions.update(retoptions)
	myHMM.produce_for_repPat(commonOptions, moreOptions)

	if not specifiedOptions["SepbamfileTemp"]==None:
		specifiedOptions["bamfile"] = (specifiedOptions["SepbamfileTemp"] % moreOptions['chr'][3:])

	if (commonOptions['SplitAndReAlign'] in [0,2]) or testall:
		start_time = time.time();
		if commonOptions['outlog'] <= M_INFO and (not specifiedOptions.has_key('thread')): print 'p2bamhmm start'
		p2bamhmm = myBAMhandler.getRepeatForGivenGene(commonOptions, specifiedOptions, moreOptions)
		memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
		if p2bamhmm==None:
			print 'ERROR None detection', moreOptions['repeatName'], moreOptions['mgloc']
			logging.error('ERROR None detection: ' + str( moreOptions['repeatName']) + ' ' + str(moreOptions['mgloc']))
		else:
			addSumForAGene(p2bamhmm, 'p2bamhmm', 2, commonOptions, specifiedOptions, moreOptions)
		end_time = time.time();
		if commonOptions['outlog'] <= M_WARNING and (not specifiedOptions.has_key('thread')): print ('p2bamhmm end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()
	if (commonOptions['SplitAndReAlign'] in [1,2]) or testall:
		start_time = time.time();
		if commonOptions['outlog'] <= M_INFO and (not specifiedOptions.has_key('thread')): print 'p2sp start'
		moreOptions['fafqfile'] = specifiedOptions["bamfile"]
		moreOptions['fafqtype'] = 'bam'
		p2sp = myRepeatReAlignment.getRepeatCounts(commonOptions, specifiedOptions, moreOptions)
		memres = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
		if p2sp==None:
			print 'ERROR None detection (sp)', moreOptions['repeatName'], moreOptions['mgloc']
			logging.error('ERROR None detection (sp): ' + str( moreOptions['repeatName']) + ' ' + str(moreOptions['mgloc']))
		else:
			addSumForAGene(p2sp, 'p2sp', 2, commonOptions, specifiedOptions, moreOptions)
		end_time = time.time();
		if commonOptions['outlog'] <= M_WARNING and (not specifiedOptions.has_key('thread')): print ('p2sp end---running time%.0f mem%d' % (end_time-start_time, memres)); sys.stdout.flush()


if __name__=="__main__":
	commonOptions = {}
	specifiedOptions = {}

	commonOptions['hg'] = 'hg38'
	commonOptions['stsBasedFolder'] = '../reference_sts/'
	commonOptions['outlog'] = M_DEBUG
	commonOptions['Patternfile'] = None
	specifiedOptions['scan_region'] = 'chr22:17551367-19651400'

	moptions = getAllMicrosatellites(commonOptions, specifiedOptions)

	mres = {}; mdetail = {}
	moreOptions = {}
	moreOptions['mres'] = mres
	moreOptions['mdetail'] = mdetail
	moreOptions['moptions'] = moptions

	micronum = filterMicrosatellites(commonOptions, specifiedOptions, moreOptions)
	allmicro = specifiedOptions['microsatellites']

	print 'Total size', micronum, allmicro[1][0].keys()


