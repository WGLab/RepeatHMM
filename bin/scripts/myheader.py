
#unique_file_id = 'myUniqueDefalt'

fastq_Sanger_quality = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

match            =  1;
mismatch         = -1;
gap_in_perf      = -1;
gap_in_read      = -10;
gap_before_after = -1



hg_reference_and_index = 'mhgversion/'

UserDefinedGenedefault = "///////"

hmm_random_rep_transit = 0.002;

template_bwamem_cmd  = 'bwa mem -k8 -W8 -r5 -A1 -B1 -O1 -E1 -L0 -t 4 %s/%s %s | samtools view -S -b | samtools sort > %s'
template_bwamem_cmd  = 'bwa mem -k8 -W8 -A1 -B1 -O1 -E1 -L0 -t 4 %s/%s %s | samtools view -S -b | samtools sort > %s'

template_bwamem_cmd  = 'bwa mem %s -A1 -B1 -O1 -E1 -L0 -t 4 %s/%s %s | samtools view -S -b | samtools sort > %s'

template_bwamem_cmd2 = 'bwa mem -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 -t 4 %s/%s %s | samtools view -S -b | samtools sort > %s'

len_isolated_repeat = 20

testall = False;
testall = True;


