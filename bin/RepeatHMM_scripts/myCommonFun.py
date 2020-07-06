
import os;
import sys;
import string;


def myWriteScanResults(specifiedOptions, mres, mdetail, procss_info='', awdefault='w'):
   scanresfolder = specifiedOptions['scanresfolder']
   if not os.path.isdir(scanresfolder):
      os.system('mkdir '+scanresfolder)
   curresfilename = scanresfolder + 'res_'+ specifiedOptions['analysis_file_id'] + procss_info+ '.txt'
   curdetailfilename = scanresfolder + 'detail_'+ specifiedOptions['analysis_file_id'] + procss_info + '.txt'
   mywrite(mres, mdetail, curresfilename, curdetailfilename, awdefault)

def myWriteScanResultsCluster(specifiedOptions, mres, mdetail, commonOptions, procss_info='', awdefault='w'):
   scanresfolder = specifiedOptions['scanresfolder']
   if not os.path.isdir(scanresfolder):
      os.system('mkdir '+scanresfolder)
   curresfilename = scanresfolder + 'res_'+ commonOptions['firsthalf_analysis_file_id'] + procss_info + commonOptions['secondhalf_analysis_file_id']
   print (curresfilename)
   curdetailfilename = scanresfolder + 'detail_'+ commonOptions['firsthalf_analysis_file_id'] + procss_info + commonOptions['secondhalf_analysis_file_id']
   mywrite(mres, mdetail, curresfilename, curdetailfilename, awdefault)

def mywrite(mres, mdetail, curresfilename, curdetailfilename, awdefault):
   msaveres = [mres, mdetail]; msavefns = [curresfilename, curdetailfilename]
   for mi in range(len(msaveres)):
      cursaveres = msaveres[mi]
      cursavefn = msavefns[mi]
      savekeys = cursaveres.keys(); savekeys.sort();
      fwriter = open(cursavefn, awdefault)
      for sk in savekeys:
        fwriter.write(str(sk))
        if isinstance(cursaveres[sk], str): #len(cursaveres[sk])==1:
           fwriter.write(' '+cursaveres[sk])
        else:
           for wi in cursaveres[sk]:
              fwriter.write(' '+str(wi))
        fwriter.write('\n')
      fwriter.close();

   if awdefault=='w':
      os.system('touch '+curresfilename+'.done')
   print ('mywrite', curresfilename+'.done')

def myReadScanResults(specifiedOptions, mres, mdetail, procss_info='', defsuf='.txt'):
   scanresfolder = specifiedOptions['scanresfolder']
   curresfilename = scanresfolder + 'res_'+ specifiedOptions['analysis_file_id'] + procss_info+ defsuf
   curdetailfilename = scanresfolder + 'detail_'+ specifiedOptions['analysis_file_id'] + procss_info + defsuf
   return myread(mres, mdetail, curresfilename, curdetailfilename)

def myReadScanResultsCluster(specifiedOptions, mres, mdetail, commonOptions, procss_info=''):
   scanresfolder = specifiedOptions['scanresfolder']
   curresfilename = scanresfolder + 'res_'+ commonOptions['firsthalf_analysis_file_id'] + procss_info + commonOptions['secondhalf_analysis_file_id']
   curdetailfilename = scanresfolder + 'detail_'+ commonOptions['firsthalf_analysis_file_id'] + procss_info + commonOptions['secondhalf_analysis_file_id']
   return myread(mres, mdetail, curresfilename, curdetailfilename)

def myread(mres, mdetail, curresfilename, curdetailfilename):
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
            print ('duplicate', curkey, cursavefn)
         cursaveres[curkey] = curv

         curline = freader.readline();
   return [mres, mdetail]
