
import os;
import sys

from scripts import findTrinucleotideRepeats


gLoc = findTrinucleotideRepeats.getDiseseGeneInRefGenomeLocation('')
glk = 'atxn3'; randtimes = 100;

if True: #False: #True:
	fqfolder = "atxn3_data/rawdata/"
	fqfile = "sam%03d.raw.fastq"
	fqstr = '_raw'
else:
	fqfolder = "atxn3_data/ccsdata/"
	fqfile = "sam%03d.ccs.fastq"
	fqstr = '_ccs'

isunsym = 0; isunsym = 1
#140         indexes = peakutils.indexes(numpy.array(y), thres=0.01) - mm
#without limiting the difference between two peaks.
for subj in range(1, 26):
	if isunsym:
		cmd = "echo 'python mySCA3_main.py -repeatgene "+glk+" --UnsymAlign 1 --align 1 --updown 18 --extend 0 --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName sca3_pcr"+str(subj)+fqstr+" --fastq "+fqfolder+(fqfile % subj)+" > logsca3/repeatgenes_"+glk+"_pcr"+str(subj)+fqstr+"_unsymalign_align_updn18_ext0_sca3.log' |  qsub -V -cwd -pe smp 1 -l h_vmem=8G -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_unsym_"+glk+'_pcr'+str(subj)
		if subj>1 and False:
			cmd += " -hold_jid trirep_unsym_"+glk+'_pcr'+str(subj-1)
	else:
		cmd = "echo 'python mySCA3_main.py -repeatgene "+glk+" --align 1 --updown 18 --extend 0 --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName sca3_pcr"+str(subj)+fqstr+" --fastq "+fqfolder+(fqfile % subj)+" > logsca3/repeatgenes_"+glk+"_pcr"+str(subj)+fqstr+"_align_updn18_ext0_sca3.log' |  qsub -V -cwd -pe smp 4 -e /home/qianliu/project/HTT_CAG_repeat/nohup -o /home/qianliu/project/HTT_CAG_repeat/nohup -N trirep_"+glk+"_pcr"+str(subj)
		if subj>1 and False:
			cmd += " -hold_jid trirep_"+glk+'_pcr'+str(subj-1)
	print cmd;
	os.system(cmd)


