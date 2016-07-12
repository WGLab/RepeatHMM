

import os;
import sys;

import trinucleotideRepeatRealSimulation

import getAlignment

if __name__=="__main__":
	ccsdatafolder = "atxn3_data/ccsdata/"
	
	if not os.path.isdir(ccsdatafolder):
		os.system("mkdir "+ccsdatafolder);

	originalcssdatafolder = 'atxn3_data/SCA3.ccs.fasta/'
	bp = getAlignment.getBasePair();
	cssids = {}

	for i in range(1, 26):
		filename = ("sam%03d.ccs.fasta" % i)
		filenameq = ("sam%03d.ccs.fastq" % i)

		fr = open(originalcssdatafolder+filename, 'r')
		fw = open(ccsdatafolder+filenameq, 'w')
		
		curline = fr.readline();
		while curline:
			curline1 = "@"+curline[1:];
			if not cssids.has_key(curline1[:-1]): cssids[curline1[:-1]] = 0
			cssids[curline1[:-1]] += 1;
			fw.write(curline1);
			
			curline2 = fr.readline();
			nalen = len(curline2)-1;
			fw.write(curline2);

			fw.write("+\n");			
			fw.write(trinucleotideRepeatRealSimulation.create_quality_for_pacbio_read(nalen, 30, 10)+'\n')
			
			
			#fw.write(curline1);
			#curline2 = getAlignment.getComplementary3(bp, curline2)

			curline = fr.readline();
		
		fr.close();
		fw.close();

	cssidskeys = cssids.keys(); cssidskeys.sort();
	for ck in cssidskeys:
		if cssids[ck]>1: print ck, cssids[ck]
	print ck
		

