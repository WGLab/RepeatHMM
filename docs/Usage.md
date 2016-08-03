# 1. Repeat count estimation for a gene

RepeatHMM is able to take a FASTQ file as input to estimate expansion count for a gene. 

## Command and Parameters:
The command and the parameters are given below:
```
usage: mySCA3_main.py [-h] [-hg HG] [--align ALIGN] [--updown UPDOWN]
                      [--extend EXTEND] [-repeatgene REPEATGENE]
                      [--UserDefinedGene USERDEFINEDGENE]
                      [--UserDefinedGeneName USERDEFINEDGENENAME]
                      [--UnsymAlign UNSYMALIGN] [--fastq FASTQ]

optional arguments:
  -h, --help            show this help message and exit
  -hg HG                The reference genome is used. Currently, only hg38 is supported
  --align ALIGN         Is unsymmetrical alignment used for error correction. 1: yes(default), 0: no
  --updown UPDOWN       Is upstream/downstream used for repeat length inference. non-0: yes(default: 18), 0: no
  --extend EXTEND       Is upstream/downstream extended as repeat region. non-0: yes, 0: no(default)
  -repeatgene REPEATGENE
                        A gene name which you want to analyze, such as HTT for huntington's disease. Default: None
  --UserDefinedGene USERDEFINEDGENE
                        The gene information defined by users. If this option is given, the default gene information will be revised. Default: ///////
  --UserDefinedGeneName USERDEFINEDGENENAME
                        The name for storing results. Default: sca3_Pcr1
  --UnsymAlign UNSYMALIGN
                        Whether unsymmetrical alignment (1) rather than bwa mem (0) would be used. Default: 1
  --fastq FASTQ         The file name for fasta sequences
```

## Example
```
For example,
        python mySCA3_main.py -repeatgene HTT;
        python mySCA3_main.py -repeatgene ATXN3;
        python mySCA3_main.py -repeatgene ATXN3 -UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1;
        python mySCA3_main.py -repeatgene all;
```
If `--UnsymAlign` is used,  `-UserDefinedGene /92070888/92072403/////` must be given and same as your PCR primer design. The two numbers are the location of the forward primer and the reverse primer.

###  with FASTQ file
```
python mySCA3_main.py -repeatgene atxn3 --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName sca3_pcr25_ccs --fastq atxn3_data/ccsdata/sam025.ccs.fastq > logsca3/repeatgenes_atxn3_pcr25_ccs_unsymalign_align_updn18_ext0_sca3.log
```

## Where to find the result:
Final results was stored in logsca3/TrinRepSca3*.log


# 2. Simulation data:

## Command and Parameters:
```
usage: trinucleotideRepeatRealSimulation_main.py [-h] [-hg HG] [--align ALIGN]
                                                 [--updown UPDOWN]
                                                 [--extend EXTEND]
                                                 [-repeatgene REPEATGENE]
                                                 [--insert_rate INSERT_RATE]
                                                 [--del_rate DEL_RATE]
                                                 [--sub_rate SUB_RATE]
                                                 [--coverage COVERAGE]
                                                 [--randTimes RANDTIMES]
                                                 [--UserDefinedGene USERDEFINEDGENE]
                                                 [--UserDefinedGeneName USERDEFINEDGENENAME]
                                                 [--UnsymAlign UNSYMALIGN]
Simulate Trinucleotide repeat for (a) known gene(s) of interests or for all genes (for hg38!!!!!!!!).

optional arguments:
  -h, --help            show this help message and exit
  -hg HG                The reference genome is used. Currently, only hg38 is supported
  --align ALIGN         Is unsymmetrical alignment used for error correction. 1: yes(Default), 0: no
  --updown UPDOWN       Is upstream/downstream used for repeat length inference. non-0: yes(Default: 18), 0: no
  --extend EXTEND       Is upstream/downstream extended as repeat region. non-0: yes, 0: no(Default)
  -repeatgene REPEATGENE
                        A gene name which you want to analyze(Default: None), such as HTT for huntington's disease.                                         'all': all known genes associated with trinucleotide repeat disorders will be analyzed;
  --insert_rate INSERT_RATE
                        Insert error rate. Default: 0.12
  --del_rate DEL_RATE   Deletion error rate. Default: 0.02
  --sub_rate SUB_RATE   Substitution error rate. Default: 0.02
  --coverage COVERAGE   The number of reads produced after mutations. Default: 300
  --randTimes RANDTIMES
                        The number of simulation times. Default: 100
  --UserDefinedGene USERDEFINEDGENE
                        The gene information defined by users. If this option is given, the default gene information will be revised. Default: ///////
  --UserDefinedGeneName USERDEFINEDGENENAME
                        The name for storing results. Default: Pcr1
  --UnsymAlign UNSYMALIGN
                        Whether unsymmetrical alignment (Default: 1) rather than bwa mem (0) would be used
```

## Example
```
For example,
        python trinucleotideRepeatRealSimulation_main.py ;
        python trinucleotideRepeatRealSimulation_main.py -repeatgene HTT;
        python trinucleotideRepeatRealSimulation_main.py -repeatgene ATXN3 ;
        python trinucleotideRepeatRealSimulation_main.py -repeatgene ATXN3 --align 1 --updown 18 --extend 0 --coverage 100 --randTimes 100 -UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1;
        python trinucleotideRepeatRealSimulation_main.py -repeatgene all;
 ```
 
 ## Where to find the result:
 Final results was stored in logsim/TrinRepSim*.log

# 3. With a BAM file
## Command and Parameters:
```
usage: findTrinucleotideRepeats_main.py [-h] [-bamfile BAMFILE] [-hg HG]
                                        [--align ALIGN] [--updown UPDOWN]
                                        [--extend EXTEND]
                                        [-repeatgene REPEATGENE]

Find trinucleotide repeats for (a) gene(s) of interests or for all genes in the test genome of the BAM file (for hg38!!!!!!!!).

optional arguments:
  -h, --help            show this help message and exit
  -bamfile BAMFILE      A BAM file storing all alignments
  -hg HG                The reference genome is used. Currently, only hg38 is supported
  --align ALIGN         Is unsymmetrical alignment used for error correction. 1: yes(Default), 0: no
  --updown UPDOWN       Is upstream/downstream used for repeat length inference. non-0: yes(Default: 18), 0: no
  --extend EXTEND       Is upstream/downstream extended as repeat region. non-0: yes, 0: no(Default )
  -repeatgene REPEATGENE
                        A gene name which you want to analyze(Default: None), such as HTT for huntington's disease.                                         'all': all known genes associated with trinucleotide repeat disorders will be analyzed;
```

## Example
```
For example,
        python findTrinucleotideRepeats_main.py -bamfile freeze4-all-merge.sort.bam;
        python findTrinucleotideRepeats_main.py -bamfile freeze4-all-merge.sort.bam -repeatgene HTT;
        python findTrinucleotideRepeats_main.py -bamfile freeze4-all-merge.sort.bam -repeatgene ATXN3 --align 1 --updown 18 --extend 0;
       python findTrinucleotideRepeats_main.py -bamfile freeze4-all-merge.sort.bam -repeatgene all;
```

 ## Where to find the result:
 Final results was stored in TrinRepDis.log




