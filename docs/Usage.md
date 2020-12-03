# The options for repeatHMM.py

There are two options below for repeatHMM.py: `BAMinput` and`FASTQinput`. Please note that users do not need to provide all parameters below. Most of them have default values.

```
usage: repeatHMM.py [-h] {BAMinput,FASTQinput,Scan} ...

Determine microsatellite repeat of interests or for all microsatellites.

positional arguments:
  {BAMinput,FASTQinput,Scan}
    BAMinput            Detect trinucleotide repeats from a BAM file
    FASTQinput          Detect trinucleotide repeats from a FASTQ file
    Scan                For scanning whole data

optional arguments:
  -h, --help            show this help message and exit

For example,
        python repeatHMM.py BAMinput: with a BAM file as input
        python repeatHMM.py FASTQinput: with a FASTQ file as input
        python repeatHMM.py Scan: for Scaning whole genome
```

## 1. Repeat count estimation for a gene

RepeatHMM is able to take a FASTQ file as input to estimate expansion count for a gene. 

### Command and Parameters:
The command and the parameters are given below:
```
usage: repeatHMM.py FASTQinput [-h] [--hg HG] [--hgfile HGFILE]
                               [--GapCorrection GAPCORRECTION]
                               [--FlankLength FLANKLENGTH]
                               [--MatchInfo MATCHINFO] [--outlog OUTLOG]
                               [--Tolerate TOLERATE] [--MinSup MINSUP]
                               [--MaxRep MAXREP] [--CompRep COMPREP]
                               [--repeatName REPEATNAME]
                               [--UserDefinedUniqID USERDEFINEDUNIQID]
                               [--Patternfile PATTERNFILE]
                               [--UserDefinedRepeat USERDEFINEDREPEAT]
                               [--SplitAndReAlign {0,1,2}]
                               [--TRFOptions TRFOPTIONS]
                               [--minTailSize MINTAILSIZE]
                               [--minRepBWTSize MINREPBWTSIZE]
                               [--RepeatTime REPEATTIME]
                               [--BWAMEMOptions BWAMEMOPTIONS]
                               [--hmm_insert_rate HMM_INSERT_RATE]
                               [--hmm_del_rate HMM_DEL_RATE]
                               [--hmm_sub_rate HMM_SUB_RATE]
                               [--SeqTech {Pacbio,Nanopore,Illumina,None}]
                               [--transitionm TRANSITIONM]
                               [--emissionm EMISSIONM] [--fastq FASTQ]

Detect trinucleotide repeats from a FASTQ file

optional arguments:
  -h, --help            show this help message and exit
  --fastq FASTQ         The file name for fasta sequences

Common options for alignment:
  --hg HG               The reference genome is used. Currently, some repeat information in hg38 and hg19 are pre-defined. The default folder is ./mhgversion/
  --hgfile HGFILE       The file name of reference genome. It could be empty and the default file name is 'hg'.fa
  --GapCorrection GAPCORRECTION
                        Is unsymmetrical alignment used for error correction of detected repeat region in reads. 1: yes(Default), 0: no
  --FlankLength FLANKLENGTH
                        Is flanking sequence used for repeat region detection. non-0: yes(Default: 30), 0: no
  --MatchInfo MATCHINFO
                        The match strategy for gap correction.

Common options:
  --outlog OUTLOG       The level for output of running information
  --Tolerate TOLERATE   Tolerate mismatch, e.g., CTG:TTG:0:2 or CTG:CTT:-2:0.
  --MinSup MINSUP       The minimum reads associated with peaks.
  --MaxRep MAXREP       The maximum repeat size. The smallest MaxRep should not be less than 14.
  --CompRep COMPREP     Whether the repeat pattern is simple ('0') or complex (nonzeo-'0': AlTlT50/C50lClT/C). For complex patterns, all patterns are required to have the same length, each position is seperated by `l` and the nucleotides at the same position is seperated by '/' where a nucleotide can by followed by a number specify relative frequency (cannot be float).

Common options for gene information:
  --repeatName REPEATNAME
                        A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. 'all': all pre-defined genes will be analyzed;
  --UserDefinedUniqID USERDEFINEDUNIQID
                        The name for storing results. Default: Pcr1. Must be different when simultaneously running RepeatHMM many times with the same setting for a same file
  --Patternfile PATTERNFILE
                        The file storing all predefined microsatellites (e.g., './reference_sts/hg38/hg38.trf.bed'). '.bed' for bed files and '.pa' for pa file. More than one files can be provided seperated by ';'.
  --UserDefinedRepeat USERDEFINEDREPEAT
                        The repeat information defined by users. If this option is given, the default repeat information will be revised. Default: ///////

Common options for re-alignment after splitting long reads using repeat regions:
  --SplitAndReAlign {0,1,2}
                        Split long reads using repeat region in long reads and re-align the non-repeat regions using BWA MEM. Default=0: not use SplitAndReAlign; 1: use SplitAndReAlign only; 2: combine 0 and 1
  --TRFOptions TRFOPTIONS
                        The options used for detecting repeat region in a read using Tandem Repeat Finder. The options are merging using _. Default='2_7_4_80_10_100_500' or '2_7_4_80_10_100'. The last parameter will be twice of the length of repeat if not given.
  --minTailSize MINTAILSIZE
                        After the split using repeat regions, discard the leftmost/rightmost non-repeat sub-sequences if they have less than --minTailSize bps. Must not be less than 10. Default=70
  --minRepBWTSize MINREPBWTSIZE
                        After the split using repeat regions, merge any two non-repeat sub-sequences if they distance is less than --minRepBWTSize bps. Must not be less than 10. Default=70
  --RepeatTime REPEATTIME
                        The minimum repeat time for a microsatellite detected by Tandem Repeat Finder. Must not be less than 2. Default=5
  --BWAMEMOptions BWAMEMOPTIONS
                        The options used for BWA MEM to align sub-reads after splitting. The options are merging using _. Default='k8_W8_r7'. For example: 'k8_W8_r7'

Common options for setting HMM matrices:
  --hmm_insert_rate HMM_INSERT_RATE
                        Insert error rate in long reads. Default: 0.12
  --hmm_del_rate HMM_DEL_RATE
                        Deletion error rate in long reads. Default: 0.02
  --hmm_sub_rate HMM_SUB_RATE
                        Substitution error rate in long reads. Default: 0.02
  --SeqTech {Pacbio,Nanopore,Illumina,None}
                        The sequencing techniques, Pacbio or Nanopore or Illumina, used to generate reads data. Default: None. Setting this option will override the setting for --hmm_insert_rate, --hmm_del_rate and --hmm_sub_rate
  --transitionm TRANSITIONM
                        User-specified transition matrix for HMM. The number of rows and columns must be the same as the 3*L+1 where L is the size of repeat unit. Probabilities in a row is separated by ',' and their sum must be 1. Probabilities of different rows are separated by ';'. Please pay more attention when providing this parameter. For CG repeat, the example of this matrix is '0.96,0.02,0,0,0,0.02,0;0,0.001,0.869,0.11,0,0,0.02;0.02,0.849,0.001,0,0.11,0.02,0;0,0.001,0.869,0.11,0,0,0.02;0.02,0.849,0.001,0,0.11,0.02,0;0.02,0.849,0.001,0,0.11,0.02,0;0.02,0.001,0.849,0.11,0,0,0.02'. Setting this option will override the setting of --hmm_insert_rate, --hmm_del_rate, --hmm_sub_rate and --SeqTech for transition matrix.
  --emissionm EMISSIONM
                        User-specified emission matrix for HMM. The number of rows must be the same as the 3*L+1 where L is the size of repeat unit and the number of columns must be 5. Probabilities in a row is separated by ',' and their sum must be 1. Probabilities of different rows are separated by ';'. Please pay more attention when providing this parameter. For CG repeat, the example of this matrix is '0.2,0.2,0.2,0.2,0.2;0.005,0.985,0.005,0.005,0;0.005,0.005,0.985,0.005,0;0.25,0.25,0.25,0.25,0;0.25,0.25,0.25,0.25,0;0.005,0.005,0.985,0.005,0;0.005,0.985,0.005,0.005,0'. Setting this option will override the setting of --hmm_insert_rate, --hmm_del_rate, --hmm_sub_rate and --SeqTech for emission matrix.
```

### Example
```
For example,
     python RepeatHMM/bin/repeatHMM.py FASTQinput --fastq XXX.fq --repeatName HTT;
     python RepeatHMM/bin/repeatHMM.py FASTQinput --fastq XXX.fq --repeatName HTT --GapCorrection 1 --FlankLength 30 ;
     python RepeatHMM/bin/repeatHMM.py FASTQinput --fastq XXX.fq --repeatName HTT --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID PCR1;
     python RepeatHMM/bin/repeatHMM.py FASTQinput --repeatName atxn3 --GapCorrection 1 --FlankLength 30 --UserDefinedUniqID sca3_pcr25_raw_test --fastq SRR5363632.fastq --outFolder logfq --Patternfile RepeatHMM/bin/reference_sts/hg38/hg38.predefined.pa --hgfile hg38/hg38.fa --hg hg38 
```
where `hg38/hg38.fa` is the hg38 human reference genome. To run the example above, users need to install `fastq-dump` and then use `fastq-dump SRR5363632`. 

### Where to find the result:
Final results was stored in logfq/RepFQ_*.log


## 2. With a BAM file
### Command and Parameters:
```
usage: repeatHMM.py BAMinput [-h] [--hg HG] [--hgfile HGFILE]
                             [--GapCorrection GAPCORRECTION]
                             [--FlankLength FLANKLENGTH]
                             [--MatchInfo MATCHINFO] [--outlog OUTLOG]
                             [--Tolerate TOLERATE] [--MinSup MINSUP]
                             [--MaxRep MAXREP] [--CompRep COMPREP]
                             [--repeatName REPEATNAME]
                             [--UserDefinedUniqID USERDEFINEDUNIQID]
                             [--Patternfile PATTERNFILE]
                             [--UserDefinedRepeat USERDEFINEDREPEAT]
                             [--SplitAndReAlign {0,1,2}]
                             [--TRFOptions TRFOPTIONS]
                             [--minTailSize MINTAILSIZE]
                             [--minRepBWTSize MINREPBWTSIZE]
                             [--RepeatTime REPEATTIME]
                             [--BWAMEMOptions BWAMEMOPTIONS]
                             [--hmm_insert_rate HMM_INSERT_RATE]
                             [--hmm_del_rate HMM_DEL_RATE]
                             [--hmm_sub_rate HMM_SUB_RATE]
                             [--SeqTech {Pacbio,Nanopore,Illumina,None}]
                             [--transitionm TRANSITIONM]
                             [--emissionm EMISSIONM]
                             [--Onebamfile ONEBAMFILE | --SepbamfileTemp SEPBAMFILETEMP]

Detect trinucleotide repeats from a BAM file

optional arguments:
  -h, --help            show this help message and exit
  --Onebamfile ONEBAMFILE
                        A BAM file storing all alignments
  --SepbamfileTemp SEPBAMFILETEMP
                        A separated BAM file template storing all alignments; separated by chromosome ids. For example, '--SepbamfileTemp'=mybam_Chr%s_sorted.bam where '%s' is chromosome id (1....22,x,y)

Common options for alignment:
  --hg HG               The reference genome is used. Currently, some repeat information in hg38 and hg19 are pre-defined. The default folder is ./mhgversion/
  --hgfile HGFILE       The file name of reference genome. It could be empty and the default file name is 'hg'.fa
  --GapCorrection GAPCORRECTION
                        Is unsymmetrical alignment used for error correction of detected repeat region in reads. 1: yes(Default), 0: no
  --FlankLength FLANKLENGTH
                        Is flanking sequence used for repeat region detection. non-0: yes(Default: 30), 0: no
  --MatchInfo MATCHINFO
                        The match strategy for gap correction.

Common options:
  --outlog OUTLOG       The level for output of running information
  --Tolerate TOLERATE   Tolerate mismatch, e.g., CTG:TTG:0:2 or CTG:CTT:-2:0.
  --MinSup MINSUP       The minimum reads associated with peaks.
  --MaxRep MAXREP       The maximum repeat size. The smallest MaxRep should not be less than 14.
  --CompRep COMPREP     Whether the repeat pattern is simple ('0') or complex (nonzeo-'0': AlTlT50/C50lClT/C). For complex patterns, all patterns are required to have the same length, each position is seperated by `l` and the nucleotides at the same position is seperated by '/' where a nucleotide can by followed by a number specify relative frequency (cannot be float).

Common options for gene information:
  --repeatName REPEATNAME
                        A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. 'all': all pre-defined genes will be analyzed;
  --UserDefinedUniqID USERDEFINEDUNIQID
                        The name for storing results. Default: Pcr1. Must be different when simultaneously running RepeatHMM many times with the same setting for a same file
  --Patternfile PATTERNFILE
                        The file storing all predefined microsatellites (e.g., './reference_sts/hg38/hg38.trf.bed'). '.bed' for bed files and '.pa' for pa file. More than one files can be provided seperated by ';'.
  --UserDefinedRepeat USERDEFINEDREPEAT
                        The repeat information defined by users. If this option is given, the default repeat information will be revised. Default: ///////

Common options for re-alignment after splitting long reads using repeat regions:
  --SplitAndReAlign {0,1,2}
                        Split long reads using repeat region in long reads and re-align the non-repeat regions using BWA MEM. Default=0: not use SplitAndReAlign; 1: use SplitAndReAlign only; 2: combine 0 and 1
  --TRFOptions TRFOPTIONS
                        The options used for detecting repeat region in a read using Tandem Repeat Finder. The options are merging using _. Default='2_7_4_80_10_100_500' or '2_7_4_80_10_100'. The last parameter will be twice of the length of repeat if not given.
  --minTailSize MINTAILSIZE
                        After the split using repeat regions, discard the leftmost/rightmost non-repeat sub-sequences if they have less than --minTailSize bps. Must not be less than 10. Default=70
  --minRepBWTSize MINREPBWTSIZE
                        After the split using repeat regions, merge any two non-repeat sub-sequences if they distance is less than --minRepBWTSize bps. Must not be less than 10. Default=70
  --RepeatTime REPEATTIME
                        The minimum repeat time for a microsatellite detected by Tandem Repeat Finder. Must not be less than 2. Default=5
  --BWAMEMOptions BWAMEMOPTIONS
                        The options used for BWA MEM to align sub-reads after splitting. The options are merging using _. Default='k8_W8_r7'. For example: 'k8_W8_r7'

Common options for setting HMM matrices:
  --hmm_insert_rate HMM_INSERT_RATE
                        Insert error rate in long reads. Default: 0.12
  --hmm_del_rate HMM_DEL_RATE
                        Deletion error rate in long reads. Default: 0.02
  --hmm_sub_rate HMM_SUB_RATE
                        Substitution error rate in long reads. Default: 0.02
  --SeqTech {Pacbio,Nanopore,Illumina,None}
                        The sequencing techniques, Pacbio or Nanopore or Illumina, used to generate reads data. Default: None. Setting this option will override the setting for --hmm_insert_rate, --hmm_del_rate and --hmm_sub_rate
  --transitionm TRANSITIONM
                        User-specified transition matrix for HMM. The number of rows and columns must be the same as the 3*L+1 where L is the size of repeat unit. Probabilities in a row is separated by ',' and their sum must be 1. Probabilities of different rows are separated by ';'. Please pay more attention when providing this parameter. For CG repeat, the example of this matrix is '0.96,0.02,0,0,0,0.02,0;0,0.001,0.869,0.11,0,0,0.02;0.02,0.849,0.001,0,0.11,0.02,0;0,0.001,0.869,0.11,0,0,0.02;0.02,0.849,0.001,0,0.11,0.02,0;0.02,0.849,0.001,0,0.11,0.02,0;0.02,0.001,0.849,0.11,0,0,0.02'. Setting this option will override the setting of --hmm_insert_rate, --hmm_del_rate, --hmm_sub_rate and --SeqTech for transition matrix.
  --emissionm EMISSIONM
                        User-specified emission matrix for HMM. The number of rows must be the same as the 3*L+1 where L is the size of repeat unit and the number of columns must be 5. Probabilities in a row is separated by ',' and their sum must be 1. Probabilities of different rows are separated by ';'. Please pay more attention when providing this parameter. For CG repeat, the example of this matrix is '0.2,0.2,0.2,0.2,0.2;0.005,0.985,0.005,0.005,0;0.005,0.005,0.985,0.005,0;0.25,0.25,0.25,0.25,0;0.25,0.25,0.25,0.25,0;0.005,0.005,0.985,0.005,0;0.005,0.985,0.005,0.005,0'. Setting this option will override the setting of --hmm_insert_rate, --hmm_del_rate, --hmm_sub_rate and --SeqTech for emission matrix.
```

### Example
```
For example,
     python repeatHMM.py BAMinput --Onebamfile XXX.bam --repeatName HTT;
     python repeatHMM.py BAMinput --Sepbamfile XXX%s.bam --repeatName HTT --GapCorrection 1 FlankLength 30;
     python repeatHMM.py BAMinput --Onebamfile freeze4-all-merge.sort.bam --repeatName HTT;
```
Please note that BAM file should be produced with the same version of 'hg' if '--hg' or '--hgfile' is specified.

### Where to find the result:
 Final results was stored in logbam/RepBAM_*.log


## 3. RepeatHMM-Scan
RepeatHMM-scan is a new module which takes the advantage of RepeatHMM to automatically scan a genome for all potential repeat regions. Together with the RepeatHMM-DB released in the assets of RepeatHMM, a genome-scale anaysis of short repeat expansions is available. 

### Command and Parameters:

```
usage: repeatHMM.py Scan [-h] [--hg HG] [--hgfile HGFILE]
                         [--GapCorrection GAPCORRECTION]
                         [--FlankLength FLANKLENGTH] [--MatchInfo MATCHINFO]
                         [--outlog OUTLOG] [--Tolerate TOLERATE]
                         [--MinSup MINSUP] [--MaxRep MAXREP]
                         [--CompRep COMPREP] [--outFolder OUTFOLDER]
                         [--repeatName REPEATNAME]
                         [--UserDefinedUniqID USERDEFINEDUNIQID]
                         [--Patternfile PATTERNFILE]
                         [--UserDefinedRepeat USERDEFINEDREPEAT]
                         [--SplitAndReAlign {0,1,2}] [--TRFOptions TRFOPTIONS]
                         [--minTailSize MINTAILSIZE]
                         [--minRepBWTSize MINREPBWTSIZE]
                         [--RepeatTime REPEATTIME]
                         [--BWAMEMOptions BWAMEMOPTIONS]
                         [--hmm_insert_rate HMM_INSERT_RATE]
                         [--hmm_del_rate HMM_DEL_RATE]
                         [--hmm_sub_rate HMM_SUB_RATE]
                         [--SeqTech {Pacbio,Nanopore,Illumina,None}]
                         [--transitionm TRANSITIONM] [--emissionm EMISSIONM]
                         [--region REGION] [--thread THREAD] [--conted CONTED]
                         [--cluster CLUSTER] [--clusterOption CLUSTEROPTION]
                         [--max_len MAX_LEN] [--StopFail STOPFAIL]
                         [--avergnum AVERGNUM] [--envset ENVSET]
                         [--repeathmmPath REPEATHMMPATH]
                         [--Onebamfile ONEBAMFILE | --SepbamfileTemp SEPBAMFILETEMP]

Scan for a whole genome or a whole region given

optional arguments:
  -h, --help            show this help message and exit
  --region REGION       The region where microsatellites are located. 'All' means all microsatellites (<10 nucleotides in repeat units), 'chrZ:X:Y' indicates a chromosome region where all microsatellites would be automatically searche. 'None': will only detect microsatellites predined in 'repeatName'
  --thread THREAD       How many additional threads are used. Default: 1
  --conted CONTED       Whether continue the running last time. Default: 0
  --cluster CLUSTER     Whether cluster computing would be used. Default: 0
  --clusterOption CLUSTEROPTION
                        Options for job submission to cluter. Default: 'qsub -V -cwd -pe smp 1 -l h_vmem=8G -q all.q,bigmem -e %s -o %s -N %s'
  --max_len MAX_LEN     The maximum length of repeat regions which will be detected. Default: 0
  --StopFail STOPFAIL   Whether stop when a job failed. Default: 0(not stop)
  --avergnum AVERGNUM   The number of repeat regions for each jobs. Default:100.
  --envset ENVSET       Whether virtual environment is created for running RepeatHMM. Default: None. Can set 'repeathmmenv'
  --repeathmmPath REPEATHMMPATH
                        The path to repeatHMM.py. Default: 'RepeatHMM/bin/'
  --Onebamfile ONEBAMFILE
                        A BAM file storing all alignments
  --SepbamfileTemp SEPBAMFILETEMP
                        A separated BAM file template storing all alignments; separated by chromosome ids. For example, '--SepbamfileTemp'=mybam_Chr%s_sorted.bam where '%s' is chromosome id (1....22,x,y)

Common options for alignment:
  --hg HG               The reference genome is used. Currently, some repeat information in hg38 and hg19 are pre-defined. The default folder is ./mhgversion/
  --hgfile HGFILE       The file name of reference genome. It could be empty and the default file name is 'hg'.fa
  --GapCorrection GAPCORRECTION
                        Is unsymmetrical alignment used for error correction of detected repeat region in reads. 1: yes(Default), 0: no
  --FlankLength FLANKLENGTH
                        Is flanking sequence used for repeat region detection. non-0: yes(Default: 30), 0: no
  --MatchInfo MATCHINFO
                        The match strategy for gap correction.

Common options:
  --outlog OUTLOG       The level for output of running information: 0:DEBUG; 1:INFO; 2:WARNING; 3:ERROR; 4:FATAL
  --Tolerate TOLERATE   Tolerate mismatch, e.g., CTG:TTG:0:2 or CTG:CTT:-2:0.
  --MinSup MINSUP       The minimum reads associated with peaks.
  --MaxRep MAXREP       The maximum repeat size. The smallest MaxRep should not be less than 14.
  --CompRep COMPREP     Whether the repeat pattern is simple ('0') or complex (nonzeo-'0': AlTlT50/C50lClT/C). For complex patterns, all patterns are required to have the same length, each position is seperated by `l` and the nucleotides at the same position is seperated by '/' where a nucleotide can by followed by a number specify relative frequency (cannot be float).
  --outFolder OUTFOLDER
                        The default folder for temporary output. Default: align/

Common options for gene information:
  --repeatName REPEATNAME
                        A gene name which you want to analyze(Default: None), such as HTT for huntington's disease. 'all': all pre-defined genes will be analyzed;
  --UserDefinedUniqID USERDEFINEDUNIQID
                        The name for storing results. Default: Pcr1. Must be different when simultaneously running RepeatHMM many times with the same setting for a same file
  --Patternfile PATTERNFILE
                        The file storing all predefined microsatellites (e.g., './reference_sts/hg38/hg38.trf.bed'). '.bed' for bed files and '.pa' for pa file. More than one files can be provided seperated by ';'.
  --UserDefinedRepeat USERDEFINEDREPEAT
                        The repeat information defined by users. If this option is given, the default repeat information will be revised. Default: ///////

Common options for re-alignment after splitting long reads using repeat regions:
  --SplitAndReAlign {0,1,2}
                        Split long reads using repeat region in long reads and re-align the non-repeat regions using BWA MEM. Default=0: not use SplitAndReAlign; 1: use SplitAndReAlign only; 2: combine 0 and 1
  --TRFOptions TRFOPTIONS
                        The options used for detecting repeat region in a read using Tandem Repeat Finder. The options are merging using _. Default='2_7_4_80_10_100_500' or '2_7_4_80_10_100'. The last parameter will be twice of the length of repeat if not given.
  --minTailSize MINTAILSIZE
                        After the split using repeat regions, discard the leftmost/rightmost non-repeat sub-sequences if they have less than --minTailSize bps. Must not be less than 10. Default=70
  --minRepBWTSize MINREPBWTSIZE
                        After the split using repeat regions, merge any two non-repeat sub-sequences if they distance is less than --minRepBWTSize bps. Must not be less than 10. Default=70
  --RepeatTime REPEATTIME
                        The minimum repeat time for a microsatellite detected by Tandem Repeat Finder. Must not be less than 2. Default=5
  --BWAMEMOptions BWAMEMOPTIONS
                        The options used for BWA MEM to align sub-reads after splitting. The options are merging using _. Default='k8_W8_r7'. For example: 'k8_W8_r7'

Common options for setting HMM matrices:
  --hmm_insert_rate HMM_INSERT_RATE
                        Insert error rate in long reads. Default: 0.12
  --hmm_del_rate HMM_DEL_RATE
                        Deletion error rate in long reads. Default: 0.02
  --hmm_sub_rate HMM_SUB_RATE
                        Substitution error rate in long reads. Default: 0.02
  --SeqTech {Pacbio,Nanopore,Illumina,None}
                        The sequencing techniques, Pacbio or Nanopore or Illumina, used to generate reads data. Default: None. Setting this option will override the setting for --hmm_insert_rate, --hmm_del_rate and --hmm_sub_rate
  --transitionm TRANSITIONM
                        User-specified transition matrix for HMM. The number of rows and columns must be the same as the 3*L+1 where L is the size of repeat unit. Probabilities in a row is separated by ',' and their sum must be 1. Probabilities of different rows are separated by ';'. Please pay more attention when providing this parameter. For CG repeat, the example of this matrix is '0.96,0.02,0,0,0,0.02,0;0,0.001,0.869,0.11,0,0,0.02;0.02,0.849,0.001,0,0.11,0.02,0;0,0.001,0.869,0.11,0,0,0.02;0.02,0.849,0.001,0,0.11,0.02,0;0.02,0.849,0.001,0,0.11,0.02,0;0.02,0.001,0.849,0.11,0,0,0.02'. Setting this option will override the setting of --hmm_insert_rate, --hmm_del_rate, --hmm_sub_rate and --SeqTech for transition matrix.
  --emissionm EMISSIONM
                        User-specified emission matrix for HMM. The number of rows must be the same as the 3*L+1 where L is the size of repeat unit and the number of columns must be 5. Probabilities in a row is separated by ',' and their sum must be 1. Probabilities of different rows are separated by ';'. Please pay more attention when providing this parameter. For CG repeat, the example of this matrix is '0.2,0.2,0.2,0.2,0.2;0.005,0.985,0.005,0.005,0;0.005,0.005,0.985,0.005,0;0.25,0.25,0.25,0.25,0;0.25,0.25,0.25,0.25,0;0.005,0.005,0.985,0.005,0;0.005,0.985,0.005,0.005,0'. Setting this option will override the setting of --hmm_insert_rate, --hmm_del_rate, --hmm_sub_rate and --SeqTech for emission matrix.

```

### Example
```
For example,
python RepeatHMM/bin/repeatHMM.py Scan --SplitAndReAlign 1 --MinSup 3 --UserDefinedUniqID chr22scan --SeqTech "Nanopore" "--Patternfile" chr22.trf.bed --cluster 1 --envset repeathmmenv --Onebamfile nanopore.minimap2.sorted.bam --hgfile GRCh38/GRCh38.fa --thread 50
```
where users need to revise 4 parameters below for their own genome: `repeathmmenv` is your repeathmm environment, `GRCh38/GRCh38.fa` is your reference genome, `nanopore.minimap2.sorted.bam` is a bam file of a nanopore data with minimap2 alignment, and `chr22.trf.bed` is potential repeat list for chr22.  `chr22.trf.bed` could be any list if it follows the format used by UCSC genome browser. For example, users can download whole-genome repeat list using `wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz`, and then after ungzip, `grep "chr22" hg38.trf.bed > chr22.trf.bed` to have `chr22.trf.bed` in this example.



