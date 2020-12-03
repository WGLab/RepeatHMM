
## Quick preparation of running RepeatHMM using conda
1. Download RepeatHMM

`git clone https://github.com/WGLab/RepeatHMM`

2. Create conda environment.
```
cd RepeatHMM
conda env create -f environment.yml
source activate repeathmmenv
cd bin/RepeatHMM_scripts/UnsymmetricPairAlignment
make
cd ../../../
```
Then, you can use run RepeatHMM by `python bin/repeatHMM.py`

If you have any error, please post them on [GitHub](https://github.com/WGLab/RepeatHMM/issues). They would also be helpful to other users.

You can also prepare the environment using the guideline below.

## Prerequisites:
	* Python 2.7
	* GCC 4.4.7
	* python packages:
		+ peakutils 1.0.3
		+ hmmlearn 0.2.1
		+ sklearn (sklearn.mixture.GaussianMixture)
		+ biopython (strongly recommend 1.66 if possible)
	  installation using pip:
	   pip install peakutils==1.0.3 hmmlearn==0.2.1 sklearn biopython==1.66
	* SWIG (see https://anaconda.org/anaconda/swig)
	* make
	* BWA MEM (see https://anaconda.org/bioconda/bwa)
	* Tandem Repeat Finder (see https://bioconda.github.io/recipes/trf/README.html)
	* samtools (see https://anaconda.org/bioconda/samtools)

## Step 1:
	* git clone https://github.com/WGLab/RepeatHMM

## Step 2:
	* go to bin/RepeatHMM_scripts/UnsymmetricPairAlignment
	* type "make" and then enter

## Usage:
 For how to use them, please refer to [Usage](https://github.com/WGLab/RepeatHMM/blob/master/docs/Usage.md)

