# The options for repeatHMM.py

```
usage: repeatHMM.py [-h] {BAMinput,FASTQinput,Simulate} ...

Determine Trinucleotide repeat for (a) gene(s) of interests or for all genes.

positional arguments:
  {BAMinput,FASTQinput,Simulate}
    BAMinput            Detect trinucleotide repeats from a BAM file
    FASTQinput          Detect trinucleotide repeats from a FASTQ file
    Simulate            For simulation

optional arguments:
  -h, --help            show this help message and exit

For example,
        python repeatHMM.py BAMinput: with a BAM file as input
        python repeatHMM.py FASTQinput: with a FASTQ file as input
        python repeatHMM.py Simulate: for simulation
```

## 1. Repeat count estimation for a gene

RepeatHMM is able to take a FASTQ file as input to estimate expansion count for a gene. 

### Command and Parameters:
The command and the parameters are given below:
```
usage: repeatHMM.py FASTQinput [-h] [--hg HG] [--hgfile HGFILE]
                               [--SeqDepth SEQDEPTH] [--RemInDel REMINDEL]
                               [--updown UPDOWN] [--extend EXTEND]
                               [--repeatgene REPEATGENE]
                               [--UserDefinedGene USERDEFINEDGENE]
                               [--UserDefinedGeneName USERDEFINEDGENENAME]
                               [--UnsymAlign UNSYMALIGN] [--fastq FASTQ]

Detect trinucleotide repeats from a FASTQ file

optional arguments:
  -h, --help            show this help message and exit
  --UnsymAlign UNSYMALIGN
                        Whether unsymmetrical alignment (1) rather than bwa mem (0) would be used. Default: 1
  --fastq FASTQ         The file name for fasta sequences

Common options for alignment:
  --hg HG               The reference genome is used. Currently, hg38 and hg19 are supported
  --hgfile HGFILE       The file name of reference genome. It could be empty and the default file name is 'hg'.fa
  --SeqDepth SEQDEPTH   The depth of the reads. For filtering in peak detection.
  --RemInDel REMINDEL   Is unsymmetrical alignment used for error correction. 1: yes(Default), 0: no
  --updown UPDOWN       Is upstream/downstream used for repeat length inference. non-0: yes(Default: 18), 0: no
  --extend EXTEND       Is upstream/downstream extended as repeat region. non-0: yes, 0: no(Default)

Common options for gene information:
  --repeatgene REPEATGENE
                        A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. 'all': all known genes associated with trinucleotide repeat disorders will be analyzed;
  --UserDefinedGene USERDEFINEDGENE
                        The gene information defined by users. If this option is given, the default gene information will be revised. Default: ///////
  --UserDefinedGeneName USERDEFINEDGENENAME
                        The name for storing results. Default: Pcr1
```

### Example
```
For example,
 python repeatHMM.py FASTQinput --fastq XXX.fq --repeatgene HTT --UnsymAlign 1;
 python repeatHMM.py FASTQinput --fastq XXX.fq --repeatgene HTT --RemInDel 1 --updown 18 --extend 0 --UnsymAlign 1;
 python repeatHMM.py FASTQinput --fastq XXX.fq --repeatgene HTT --RemInDel 1 --updown 18 --extend 0 --UnsymAlign 1 --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1;
```
If `--UnsymAlign` is used,  `-UserDefinedGene /92070888/92072403/////` must be given and same as your PCR primer design. The two numbers are the location of the forward primer and the reverse primer.

####  with FASTQ file
```
 python repeatHMM.py FASTQinput --repeatgene atxn3 --UnsymAlign 1 --RemInDel 1 --updown 18 --extend 0 --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName sca3_pcr25_raw_test --fastq atxn3_data/rawdata/sam025.raw.fastq
```

### Where to find the result:
Final results was stored in logfq/TrinRepFQ*.log


## 2. Simulation data:

### Command and Parameters:
```
usage: repeatHMM.py Simulate [-h] [--hg HG] [--hgfile HGFILE]
                             [--SeqDepth SEQDEPTH] [--RemInDel REMINDEL]
                             [--updown UPDOWN] [--extend EXTEND]
                             [--repeatgene REPEATGENE]
                             [--UserDefinedGene USERDEFINEDGENE]
                             [--UserDefinedGeneName USERDEFINEDGENENAME]
                             [--insert_rate INSERT_RATE] [--del_rate DEL_RATE]
                             [--sub_rate SUB_RATE] [--coverage COVERAGE]
                             [--randTimes RANDTIMES] [--UnsymAlign UNSYMALIGN]

Simulate for a given gene

optional arguments:
  -h, --help            show this help message and exit
  --insert_rate INSERT_RATE
                        Insert error rate. Default: 0.12
  --del_rate DEL_RATE   Deletion error rate. Default: 0.02
  --sub_rate SUB_RATE   Substitution error rate. Default: 0.02
  --coverage COVERAGE   The number of reads produced after mutations. Default: 300
  --randTimes RANDTIMES
                        The number of simulation times. Default: 100
  --UnsymAlign UNSYMALIGN
                        Whether unsymmetrical alignment (1) rather than bwa mem (0) would be used. Default: 1

Common options for alignment:
  --hg HG               The reference genome is used. Currently, hg38 and hg19 are supported
  --hgfile HGFILE       The file name of reference genome. It could be empty and the default file name is 'hg'.fa
  --SeqDepth SEQDEPTH   The depth of the reads. For filtering in peak detection.
  --RemInDel REMINDEL   Is unsymmetrical alignment used for error correction. 1: yes(Default), 0: no
  --updown UPDOWN       Is upstream/downstream used for repeat length inference. non-0: yes(Default: 18), 0: no
  --extend EXTEND       Is upstream/downstream extended as repeat region. non-0: yes, 0: no(Default)

Common options for gene information:
  --repeatgene REPEATGENE
                        A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. 'all': all known genes associated with trinucleotide repeat disorders will be analyzed;
  --UserDefinedGene USERDEFINEDGENE
                        The gene information defined by users. If this option is given, the default gene information will be revised. Default: ///////
  --UserDefinedGeneName USERDEFINEDGENENAME
                        The name for storing results. Default: Pcr1
```

### Example
```
For example,
 python repeatHMM.py Simulate --repeatgene HTT --UnsymAlign 1;
 python repeatHMM.py Simulate --repeatgene HTT --RemInDel 1 --updown 18 --extend 0 --UnsymAlign 1;
 python repeatHMM.py Simulate --repeatgene HTT --RemInDel 1 --updown 18 --extend 0 --UnsymAlign 1 --coverage 100 --randTimes 100 -UserDefinedGene /92070888/92072403///// --UserDefinedGeneName PCR1;
 python repeatHMM.py Simulate --repeatgene atxn3 --UnsymAlign 1 --RemInDel 1 --updown 18 --extend 0 --randTimes 100 --coverage 50
 python repeatHMM.py Simulate --repeatgene atxn3 --UnsymAlign 1 --RemInDel 1 --updown 18 --extend 0 --randTimes 100 --coverage 50 --UserDefinedGene /92070888/92072403///// --UserDefinedGeneName pcr1
 ```
 
 ### Where to find the result:
 Final results was stored in logsim/TrinRepSim*.log

## 3. With a BAM file
### Command and Parameters:
```
usage: repeatHMM.py BAMinput [-h] [--hg HG] [--hgfile HGFILE]
                             [--SeqDepth SEQDEPTH] [--RemInDel REMINDEL]
                             [--updown UPDOWN] [--extend EXTEND]
                             [--repeatgene REPEATGENE]
                             [--UserDefinedGene USERDEFINEDGENE]
                             [--UserDefinedGeneName USERDEFINEDGENENAME]
                             [--Onebamfile ONEBAMFILE | --SepbamfileTemp SEPBAMFILETEMP]

Detect trinucleotide repeats from a BAM file

optional arguments:
  -h, --help            show this help message and exit
  --Onebamfile ONEBAMFILE
                        A BAM file storing all alignments
  --SepbamfileTemp SEPBAMFILETEMP
                        A separated BAM file template storing all alignments; separated by chromosome ids. For example, '--SepbamfileTemp'=mybam_Chr%s_sorted.bam where '%s' is chromosome id (1....22,x,y)

Common options for alignment:
  --hg HG               The reference genome is used. Currently, hg38 and hg19 are supported
  --hgfile HGFILE       The file name of reference genome. It could be empty and the default file name is 'hg'.fa
  --SeqDepth SEQDEPTH   The depth of the reads. For filtering in peak detection.
  --RemInDel REMINDEL   Is unsymmetrical alignment used for error correction. 1: yes(Default), 0: no
  --updown UPDOWN       Is upstream/downstream used for repeat length inference. non-0: yes(Default: 18), 0: no
  --extend EXTEND       Is upstream/downstream extended as repeat region. non-0: yes, 0: no(Default)

Common options for gene information:
  --repeatgene REPEATGENE
                        A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. 'all': all known genes associated with trinucleotide repeat disorders will be analyzed;
  --UserDefinedGene USERDEFINEDGENE
                        The gene information defined by users. If this option is given, the default gene information will be revised. Default: ///////
  --UserDefinedGeneName USERDEFINEDGENENAME
                        The name for storing results. Default: Pcr1
```

### Example
```
For example,
 python repeatHMM.py BAMinput --bamfile XXX.bam --repeatgene HTT;
 python repeatHMM.py BAMinput --bamfile XXX.bam --repeatgene HTT --RemInDel 1 --updown 18 --extend 0;
 python repeatHMM.py BAMinput --bamfile freeze4-all-merge.sort.bam --repeatgene HTT;
```
Please note that BAM file should be produced with the same version of 'hg' if '--hg' or '--hgfile' is specified.

### Where to find the result:
 Final results was stored in logbam/TrinRepBAM*.log




