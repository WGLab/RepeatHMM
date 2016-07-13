
import os;
import sys

from scripts import findTrinucleotideRepeats


gLoc = findTrinucleotideRepeats.getDiseseGeneInRefGenomeLocation('')

glk = 'atxn3'; randtimes = 100;

isunsym = False; isunsym = 1
if isunsym:
	ispcr = True; ispcr = False;
	if ispcr:
		allcoverages = range(20, 100, 20);
		allcoverages.extend(range(100, 1000, 100));
		allcoverages.extend(range(1000, 10001, 300))
	else:
		allcoverages = range(20, 100, 20);
		allcoverages.extend(range(100, 1000, 100));
		#allcoverages.extend(range(1000, 5001, 300))
		allcoverages.extend(range(1000, 10001, 300))
	print allcoverages

	for coverage in allcoverages:
                if (not ispcr):
                        cmd = "echo 'python trinucleotideRepeatRealSimulation_main.py -repeatgene "+glk+" --UnsymAlign 1 --align 1 --updown 18 --extend 0 --randTimes "+str(randtimes)+" --coverage "+str(coverage)+" > logsim/repeatgenes_"+glk+"_unsymalign_align_updn18_nonext_errorDef_coverage"+str(coverage)+"_time"+str(randtimes)+".log' | qsub -V -cwd -pe smp 1 -l h_vmem=7G -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_unsym_"+glk+"_"+str(coverage)+" -hold_jid " + "trirep_unsym_"+glk+"_"+str(coverage)+'_pcr1'+" -hold_jid trirep_"+glk+"_"+str(coverage)+"_align_updn120_ext6_errorDef_pcr1"
                else:
                        cmd = "echo 'python trinucleotideRepeatRealSimulation_main.py -repeatgene "+glk+" --UnsymAlign 1 --align 1 --updown 18 --extend 0 --randTimes "+str(randtimes)+" --coverage "+str(coverage)+" --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName pcr1 > logsim/repeatgenes_"+glk+"_pcr1_unsymalign_align_updn18_nonext_errorDef_coverage"+str(coverage)+"_time"+str(randtimes)+".log' |  qsub -V -cwd -pe smp 1 -l h_vmem=7G -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_unsym_"+glk+"_"+str(coverage)+'_pcr1'
			#on Jul 7, 2016, only work when coverage < 7000
                print cmd;

                os.system(cmd)
	sys.exit(0)

allcoverages = range(20, 100, 20);
allcoverages.extend(range(100, 1000, 100));
allcoverages.extend(range(1000, 10001, 300))
print allcoverages

istest = False;
if istest:
	allcoverages = range(4900, 4950, 100);

isalign = 1; ispcr = True;

for coverage in allcoverages:
	if isalign:
                if (not ispcr):
                        cmd = "echo 'python trinucleotideRepeatRealSimulation_main.py -repeatgene "+glk+" --align 1 --updown 120 --extend 6 --randTimes "+str(randtimes)+" --coverage "+str(coverage)+" > logsim/repeatgenes_"+glk+"_align_updn120_ext6_errorDef_coverage"+str(coverage)+"_time"+str(randtimes)+".log' | qsub -V -cwd -pe smp 4 -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_"+glk+"_"+str(coverage)+"_align_updn120_ext6_errorDef"
                else:
                        cmd = "echo 'python trinucleotideRepeatRealSimulation_main.py -repeatgene "+glk+" --align 1 --updown 120 --extend 6 --randTimes "+str(randtimes)+" --coverage "+str(coverage)+" --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName pcr1 > logsim/repeatgenes_"+glk+"_pcr1_align_updn120_ext6_errorDef_coverage"+str(coverage)+"_time"+str(randtimes)+".log' |  qsub -V -cwd -pe smp 4 -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_"+glk+"_"+str(coverage)+"_align_updn120_ext6_errorDef_pcr1 -hold_jid trirep_"+glk+"_"+str(coverage)+"_align_updn120_ext6_errorDef"
                print cmd;

                os.system(cmd)
	elif (not istest):
		if (not ispcr):
			cmd = "echo 'python trinucleotideRepeatRealSimulation_main.py -repeatgene "+glk+" --align 0 --updown 0 --extend 6 --randTimes "+str(randtimes)+" --coverage "+str(coverage)+" > logsim/repeatgenes_"+glk+"_non_align_updn0_ext6_errorDef_coverage"+str(coverage)+"_time"+str(randtimes)+".log' | qsub -V -cwd -pe smp 4 -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_"+glk+"_"+str(coverage)
		else:
			cmd = "echo 'python trinucleotideRepeatRealSimulation_main.py -repeatgene "+glk+" --align 0 --updown 0 --extend 6 --randTimes "+str(randtimes)+" --coverage "+str(coverage)+" --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName pcr1 > logsim/repeatgenes_"+glk+"_pcr1_non_align_updn0_ext6_errorDef_coverage"+str(coverage)+"_time"+str(randtimes)+".log' |  qsub -V -cwd -pe smp 4 -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_"+glk+"_"+str(coverage)+'_pcr1'
		print cmd;
		
		os.system(cmd)
		
	else:
		cmd = "echo 'python trinucleotideRepeatRealSimulation_main.py -repeatgene "+glk+" --align 0 --updown 0 --extend 3 --randTimes "+str(randtimes)+" --coverage "+str(coverage)+" > logsim/repeatgenes_"+glk+"_non_align_updn0_ext3_errorDef_coverage"+str(coverage)+"_time"+str(randtimes)+".log' | qsub -V -cwd -pe smp 4 -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_"+glk+"_"+str(coverage)
		print cmd;
		
		cmd = "echo 'python trinucleotideRepeatRealSimulation_main.py -repeatgene "+glk+" --align 0 --updown 0 --extend 3 --randTimes "+str(randtimes)+" --coverage "+str(coverage)+" --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1 > logsim/repeatgenes_"+glk+"_pcr1_non_align_updn0_ext3_errorDef_coverage"+str(coverage)+"_time"+str(randtimes)+".log' |  qsub -V -cwd -pe smp 4 -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_"+glk+"_"+str(coverage)+'_pcr1'
		print cmd;



