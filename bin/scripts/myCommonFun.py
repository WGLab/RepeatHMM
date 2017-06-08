
import os;
import sys;
import string;


def myWriteScanResults(specifiedOptions, mres, mdetail, procss_info='', awdefault='w'):
   scanresfolder = specifiedOptions['scanresfolder']
   if not os.path.isdir(scanresfolder):
      os.system('mkdir '+scanresfolder)
   curresfilename = scanresfolder + 'res_'+ specifiedOptions['analysis_file_id'] + procss_info+ '.txt'
   curdetailfilename = scanresfolder + 'detail_'+ specifiedOptions['analysis_file_id'] + procss_info + '.txt'
   msaveres = [mres, mdetail]; msavefns = [curresfilename, curdetailfilename]
   for mi in range(len(msaveres)):
      cursaveres = msaveres[mi]
      cursavefn = msavefns[mi]
      savekeys = cursaveres.keys(); savekeys.sort();
      fwriter = open(cursavefn, awdefault)
      for sk in savekeys:
        fwriter.write(str(sk))
        for wi in cursaveres[sk]:
           fwriter.write(' '+str(wi))
        fwriter.write('\n')
      fwriter.close();

   if awdefault=='w':
      os.system('touch '+curresfilename+'.done')

def myReadScanResults(specifiedOptions, mres, mdetail, procss_info='', defsuf='.txt'):
	scanresfolder = specifiedOptions['scanresfolder']
	curresfilename = scanresfolder + 'res_'+ specifiedOptions['analysis_file_id'] + procss_info+ defsuf
	curdetailfilename = scanresfolder + 'detail_'+ specifiedOptions['analysis_file_id'] + procss_info + defsuf
	msaveres = [mres, mdetail]; msavefns = [curresfilename, curdetailfilename]
	for mi in range(len(msaveres)):
		cursaveres = msaveres[mi]
		cursavefn = msavefns[mi]
		freader = open(cursavefn, 'r')
		
		curline = freader.readline();
		while curline:
			curline = string.strip(curline);
			space1ind = curline.index(' ');
			curkey = curline[:space1ind]
			curv = curline[space1ind+1:]

			if cursaveres.has_key(curkey):
				print 'duplicate', curkey, cursavefn
			cursaveres[curkey] = curv

			curline = freader.readline();
	return [mres, mdetail]
