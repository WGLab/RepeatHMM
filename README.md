# RepeatHMM: estimation of repeat counts on microsatellites from long-read sequencing data

RepeatHMM is a novel computational tool to detect any microsatellites (including trinucleotide repeats in trinucleotide repeat disorders (TRD)) from given long reads for a subject of interests. It is able to accurately estimate estimate expansion counts according to the evaluation performance on both simulation data and real data. It is user friendly and easy to install and use.

## Features

* Accurate and efficient estimation of repeat counts from long-read sequencing data

* Analysis of all types of simple repeats

* Prefined models are included for more than 10 well known trinucleotide repeats: AFF2, AR, ATN1, ATXN1, ATXN2, ATXN3, ATXN7, ATXN8OS, CACNA1A, DMPK, FMR1, FXN, HTT, PPP2R2B, TBP

* Easy to install and use

## Methodology of RepeatHMM

RepeatHMM takes a set of reads as input, uses a split-and-align strategy to improve alignments, performs error correction, and leverages a hidden Markov model (HMM) and a peak calling algorithm based on Gaussian mixture model to infer repeat counts. RepeatHMM allows users to specify error parameters of the sequencing experiments, thus automatically producing transition and emission matrices for HMM and allowing the analysis of both PacBio and Oxford Nanopore data. 

RepeatHMM was evaluated on both random simulation and PCR-based simulation for long reads containing CAG repeats, and also on real datasets of ATXN3 for SCA3 of ATXN10  for SCA10. The results demonstrated that our tool was able to accurately estimate expansion counts from long reads.

## Inputs of RepeatHMM

RepeatHMM takes long reads from a subject as input, and can also take a BAM file as input to find more than 10 predefined trinucleotide repeats or a gene given by users, after all reads were well aligned to a reference genome. 

## Usage

Please refer to [Usage](https://github.com/WGLab/RepeatHMM/blob/master/docs/Usage.md) for how to use RepeatHMM.

## Revision History

For release history, please visit [here](https://github.com/WGLab/RepeatHMM/releases). For details, please go [here](https://github.com/WGLab/RepeatHMM/blob/master/README.md).

## Contact

If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/WGLab/RepeatHMM/issues). They would also be helpful to other users. 

## Reference

Liu Q, Zhang P, Wang D, Gu W, Wang K. Interrogating the "unsequenceable" genomic trinucleotide repeat disorders by long-read sequencing. Genome Medicine, in press, 2017
