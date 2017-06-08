
#fastq_Sanger_quality = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

M_DEBUG = 0;
M_INFO = 1;
M_WARNING = 2;
M_ERROR = 3;
M_FATAL = 4;

cur_M_STAT = M_WARNING;

FATAL_key = 'FATAL'

#match            =  3 # 5 #
#mismatch         = -2;
#gap_in_perf      = -2;
#gap_in_read      = -15;
#gap_before_after = -1
##no toerlate
##5, -15;  5, -11; 3, -11; 3, -15;   
## pattern tolerate 3, -15;
## pattern tolerate2 3, -15;

min_flank_len = 10;

hg_reference_and_index = 'mhgversion/'
stsBasedFolder = 'reference_sts/'
UserDefinedRepeatdefault = "///////"

hmm_random_rep_transit = 0.002;

mthreads = '4';
mthreads = '1'

template_bwamem_cmd  = 'bwa mem %s -A1 -B1 -O1 -E1 -L0 -t '+mthreads+' -v 2 %s/%s %s | samtools view -S -b | samtools sort > %s'
template_bwamem_cmd2 = 'bwa mem -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 -t '+mthreads+' -v 2 %s/%s %s | samtools view -S -b | samtools sort > %s'


len_isolated_repeat = 20

testall = False;
testall = True;


