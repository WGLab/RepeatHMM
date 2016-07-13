

import os;
import sys;

import string;

from scripts import trinucleotideRepeatRealSimulation

def getMean_N50(lenlist, totbase):
	lenlist.sort();
	
	llen = len(lenlist);
	mid = llen/2
	if llen%2==1: mean = lenlist[mid]
	else:
		mean = (lenlist[mid] + lenlist[mid-1])/2

	bashalf = totbase/2
	curhalf = 0; rev_ind = -1;
	while curhalf<bashalf and len(lenlist)+rev_ind>0:
		curhalf += lenlist[rev_ind]
		if curhalf>=bashalf: break;
		else: rev_ind -= 1;
	
	if curhalf>bashalf: n50 = lenlist[rev_ind]
	elif curhalf==bashalf: n50 = (lenlist[rev_ind]+lenlist[rev_ind-1])/2
	else: print 'Error could not find n50', len(lenlist), rev_ind

	return [n50, mean, len(lenlist)+rev_ind];


if __name__=="__main__":
	ccsdatafolder = "atxn3_data/rawdata/"
	
	if not os.path.isdir(ccsdatafolder):
		os.system("mkdir "+ccsdatafolder);

	originalcssdatafolder = 'atxn3_data/SCA3.ccs.fasta/'
	cssids = {}; fws = []; errstr = ''; lessthan1500 = []; lessthan1400=[]; lessthan1300=[];
	mfiecount = []; totbase=0; toalenlist = []; lenlistall = []
	for i in range(1, 26):
		lessthan1500.append(0); lessthan1400.append(0); lessthan1300.append(0);
		lenlistall.append([])
		mfiecount.append([0,0, 0, 0]); lenlist1 = [];
		filename = ("sam%03d.ccs.fasta" % i)
		fr = open(originalcssdatafolder+filename, 'r')
		curline = fr.readline();
		while curline:
			mfiecount[-1][0] += 1
			
			lsp = curline.split('/')
			if lsp[0][0]=='>' or lsp[0][0]=='<' or lsp[0][0]=='@':
				curkey = lsp[0][1:]+'/'+lsp[1]
			else:
				curkey = lsp[0]+'/'+lsp[1]
				print 'Warning', curline, lsp, curkey
			if not cssids.has_key(curkey): cssids[curkey] = []
			cssids[curkey].append(i);
			if len(cssids[curkey])>1:
				errstr +=  'Error, more file with the same id: ' + curline
				for fid in cssids[curkey]: errstr += ' '+str(fid);
				errstr += '\n'

			curline = string.strip(fr.readline());

			lenlist1.append(len(curline))
			mfiecount[-1][2] += lenlist1[-1]
			if lenlist1[-1]<1500:
				print i, lenlist1[-1]
				lessthan1500[-1] += 1
				if lenlist1[-1]<1497:
					lessthan1400[-1] += 1
					if lenlist1[-1]<1494: lessthan1300[-1] += 1
			totbase += lenlist1[-1]
			toalenlist.append(lenlist1[-1])

			curline = fr.readline();
		mfiecount[-1].append(getMean_N50(lenlist1, mfiecount[-1][2]));
		fr.close();
		filenameq = ("sam%03d.raw.fastq" % i)
		fws.append(open(ccsdatafolder+filenameq, 'w'));
	if not errstr=='':
		print errstr
		for i in range(len(fws)): fws[i].close();
		os.exit(1);
	
	for lt1300_ind in range(len(lessthan1300)):
		if lessthan1300[lt1300_ind]>0 or lessthan1500[lt1300_ind]>0 or lessthan1400[lt1300_ind]>0: 
			print 'lessthan1300 in CCS(1494, 1497, 1500): ', " {:2,}".format(lt1300_ind+1), ':', " {:5,}".format(lessthan1300[lt1300_ind]), " {:5,}".format(lessthan1400[lt1300_ind]), " {:5,}".format(lessthan1500[lt1300_ind])
			#if lessthan1500[lt1300_ind]==0 or lessthan1400[lt1300_ind]==0: print 'Warning in CSS', lt1300_ind+1, lessthan1300[lt1300_ind], lessthan1400[lt1300_ind],lessthan1500[lt1300_ind]

		lessthan1500[lt1300_ind] = 0; lessthan1400[lt1300_ind] = 0; lessthan1300[lt1300_ind] = 0;
		
	
	mfiecount.append([len(toalenlist), 0, totbase, 0])
	mfiecount[-1].append(getMean_N50(toalenlist, totbase))

	totbase=0; toalenlist = []; 
	cssidskeys = cssids.keys(); cssidskeys.sort();
	allfr = open('atxn3_data/all.fq', 'r')
	acurl = allfr.readline(); nokeynum = 0; total = 0;
	while acurl:
		total += 1; 
		lsp = acurl.split('/')
		curkey = lsp[0][1:]+'/'+lsp[1]
		if cssids.has_key(curkey):
			curfw = fws[cssids[curkey][0]-1]
			mfiecount[cssids[curkey][0]-1][1] += 1
			#@
			curfw.write(acurl);
			#na
			acurl = allfr.readline();

			lenlistall[cssids[curkey][0]-1].append(len(string.strip(acurl)))
			mfiecount[cssids[curkey][0]-1][3] += lenlistall[cssids[curkey][0]-1][-1]
			if lenlistall[cssids[curkey][0]-1][-1]<1500:
				lessthan1500[cssids[curkey][0]-1] += 1
				if lenlistall[cssids[curkey][0]-1][-1]<1400:
					lessthan1400[cssids[curkey][0]-1] += 1
					if lenlistall[cssids[curkey][0]-1][-1]<1300: lessthan1300[cssids[curkey][0]-1] += 1
			totbase += lenlistall[cssids[curkey][0]-1][-1]
			toalenlist.append(lenlistall[cssids[curkey][0]-1][-1])

			curfw.write(acurl);
			#+
			acurl = allfr.readline();
			curfw.write(acurl);
			#quality
			acurl = allfr.readline();
			curfw.write(acurl);
			
			acurl = allfr.readline();
		else:
			#print 'Error no key for ', acurl[:-1], '\t', curkey, '\t', cssidskeys[0]
			acurl = allfr.readline();
			acurl = allfr.readline();
			acurl = allfr.readline();
			acurl = allfr.readline();
			nokeynum += 1;

	for mfc in range(len(lenlistall)):
		mfiecount[mfc].append(getMean_N50(lenlistall[mfc], mfiecount[mfc][3]));

	mfiecount[-1][1] = len(toalenlist); mfiecount[-1][3] = totbase
	mfiecount[-1].append(getMean_N50(toalenlist, totbase))
	
	allfr.close();
	for i in range(len(fws)):
		fws[i].close();

	for lt1300_ind in range(len(lessthan1300)):
		if lessthan1300[lt1300_ind]>0 or lessthan1500[lt1300_ind]>0 or lessthan1400[lt1300_ind]>0: 
			print 'lessthan1300 in raw(1300, 1400, 1500):', " {:2,}".format(lt1300_ind+1), ':', " {:5,}".format(lessthan1300[lt1300_ind]), " {:5,}".format(lessthan1400[lt1300_ind]), " {:5,}".format(lessthan1500[lt1300_ind])
			#if lessthan1500[lt1300_ind]==0 or lessthan1400[lt1300_ind]==0: print 'Warning in CSS', lt1300_ind+1, lessthan1300[lt1300_ind], lessthan1400[lt1300_ind],lessthan1500[lt1300_ind]

		lessthan1500[lt1300_ind] = 0; lessthan1400[lt1300_ind] = 0; lessthan1300[lt1300_ind] = 0;	

	print 'nokeynum', (" {:,}".format(nokeynum)), 'total', (" {:,}".format(total))
	for mfc in range(len(mfiecount)):
		print ('%2d: [' % (mfc+1)),
		for it in mfiecount[mfc]:
			if isinstance(it, list):
				print ' [',
				for dei in it:
					print (" {:6,}".format(dei)), 	
				print ' ]',
			else:
				print (" {:6,}".format(it)),
		print ' ]'
		if (mfc+1)%5==0: print ''


	print '\nCCS'
	for mfc in range(len(mfiecount)):
		print ('%2d:' % (mfc+1)),
		print (" {:6,}".format(mfiecount[mfc][0])),
		print (" {:6,}".format(mfiecount[mfc][2])),
		print (" {:6,}".format(mfiecount[mfc][4][0])),
		print (" {:6,}".format(mfiecount[mfc][4][1]))
	print '\nraw'
	for mfc in range(len(mfiecount)):
		print ('%2d:' % (mfc+1)),
		print (" {:6,}".format(mfiecount[mfc][1])),
		print (" {:6,}".format(mfiecount[mfc][3])),
		print (" {:6,}".format(mfiecount[mfc][5][0])),
		print (" {:6,}".format(mfiecount[mfc][5][1]))
	
