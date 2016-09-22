
import os;
import sys;
import math;
import numpy as np

from scripts import trinucleotideRepeatRealSimulation

if __name__=='__main__':
	glk = 'atxn3'; randtimes = 100;

	allcoverages = range(20, 100, 20);
	allcoverages.extend(range(100, 1000, 100));
	#allcoverages.extend(range(1000, 10001, 300))
	
	allpcrstr = ['_defpcr1', ''];
	
	for pcrstr in allpcrstr:
		allres = {};

		if pcrstr=='':
			allcoverages = range(10, 100, 10);
			allcoverages.extend(range(100,  1000, 100));
			#allcoverages.extend(range(1000, 10001, 300))
			allcoverages.extend(range(1000, (3000+1), 1000));
		else:
			allcoverages = range(10, 100, 10);
			allcoverages.extend(range(100,  1000, 100));
			#allcoverages.extend(range(1000, 10001, 300))
			allcoverages.extend(range(1000, (3000+1), 1000));
	
		for coverage in allcoverages:
			#curf = "./sim_res/atxn3"+pcrstr+"_ins0.12_del0.02_sub0.02_cov"+str(coverage)+"_unsymalign_align_updn18_non_ext_times"+str(randtimes)+".new.txt"
			curf = "./sim_res/atxn3"+pcrstr+"_ins0.12_del0.02_sub0.02_cov"+str(coverage)+"_unsymalign_align_updn18_non_ext_times"+str(randtimes)+".txt"


			#print curf;
			#if (pcrstr=='') and coverage==200: print curf
			allres[coverage] = []
			allpredsim = np.array(trinucleotideRepeatRealSimulation.readList(curf));
			for i in range(len(allpredsim)):
				dif1 = math.fabs(allpredsim[i][0] - allpredsim[i][2])
				dif2 = math.fabs(allpredsim[i][1] - allpredsim[i][3])

				difA = dif2
				if dif1<dif2: difA = dif1
			
				dif3 = math.fabs(allpredsim[i][0] - allpredsim[i][3])
				dif4 = math.fabs(allpredsim[i][1] - allpredsim[i][2])
			
				difB = dif3
				if dif3>dif4: difB = dif4
			
				if difB+5 < difA: 
					#pass 
					allpredsim[i][2], allpredsim[i][3] = allpredsim[i][3], allpredsim[i][2]
			
				if (dif1<4 and dif2<4) or (dif3<4 and dif4<4): pass
				else: 
					#if (pcrstr=='') and coverage==200:
					#	print '\t', coverage, allpredsim[i]
					#else:
						pass # print '\t', coverage, allpredsim[i]
				if (pcrstr=='') and coverage==200:
					pass #print coverage, allpredsim[i]

		
			allpredsim = allpredsim.T		
	
			rmse1 = np.sqrt(np.mean((allpredsim[0,]-allpredsim[2,])**2))
			nozero1 = np.sum(allpredsim[0,]-allpredsim[2,]!=0)
			rmse2 = np.sqrt(np.mean((allpredsim[1,]-allpredsim[3,])**2))
			nozero2 = np.sum(allpredsim[1,]-allpredsim[3,]!=0)
	
			rmse1b = np.sqrt(np.mean((allpredsim[0,]-allpredsim[4,])**2))
			rmse2b = np.sqrt(np.mean((allpredsim[1,]-allpredsim[5,])**2))
	
			allres[coverage].append(rmse1)
			allres[coverage].append(rmse2)
			allres[coverage].append(rmse1b); allres[coverage].append(rmse2b)
			if (pcrstr=='') and coverage==200: print allres[coverage], '\n'

		for coverage in allcoverages:
			#print ("%6d: %7.3f, %7.3f" % (coverage, allres[coverage][0], allres[coverage][1]))
			print ("%d %.3f %.3f %.3f %.3f" % (coverage, allres[coverage][0], allres[coverage][1], allres[coverage][2], allres[coverage][3]))
		print ''
