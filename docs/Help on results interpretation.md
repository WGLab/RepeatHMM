# Some help document to run `repeatHMM.py`

## Where to find the result and how to interpret it

`repeatHMM.py FASTQinput` will store the result in `logfq/RepFQ_*.log`, while `repeatHMM.py BAMinput` stores result in `logbam/RepBAM_*.log`. 

Usually, you will have some line in the result file like https://github.com/WGLab/RepeatHMM/issues/2:
```
p2bamhmm ['atxn3', 14.0, [27, 27], 'allocr:11:1, 14:1, 15:1, 16:1, 18:1, 19:4, 20:195, 21:6, 22:3, 23:2, 24:7, 25:76, 26:626, 27:519, 28:51, 29:7, 30:3, 31:3, 32:1, 34:1, 35:2, 36:1, 37:2, 39:2, 40:3, 41:2, 42:3, 43:7, 44:1, 45:4, 46:6, 47:3, 48:6, 49:2, 50:4, 51:3, 52:4, 53:1, 54:2, 57:1, 62:1, 63:1, 64:1, 67:2, 69:1, 70:2, 71:2, 72:2,', 1580, 205]
p2sp ['atxn3', 14.0, [27, 69], 'allocr:20:1, 22:1, 23:2, 24:7, 25:76, 26:625, 27:520, 28:51, 29:6, 30:3, 31:3, 32:1, 35:2, 36:1, 37:2, 39:2, 40:3, 41:2, 42:3, 43:5, 45:4, 46:6, 47:4, 48:5, 49:3, 50:6, 51:4, 52:6, 53:2, 54:6, 55:4, 56:1, 57:2, 58:4, 59:4, 60:4, 61:3, 62:2, 63:3, 64:2, 65:2, 66:5, 67:11, 68:24, 69:35, 70:43, 71:35, 72:27, 73:7, 74:4, 75:5, 78:1, 83:1,', 1591, 2]
```
p2sp's results are recommended which usually are better than p2bamhmm's results. p2bamhmm is a simple version of HMM and it is good for some situations where one need short running time but pay less attention to the accuracy. [XX, YY] is the estimated allele size, followed by a distribution of all estimated repeat sizes.

## `INFO: Warning!!!! could not find optimized model`
Sometimes, you might find the warning above in the result file. Do not be panic! This is because Gaussian Mixture model is used in RepeatHMM to infer repeat count. Thus, enough reads are needed to support the estimation. When there are no enough supporting reads and each estimated repeat count is only supported by very few reads. Gaussian Mixture model thus cannot find a better model to fit the data and give the warning. Thus, the suitable coverage is necessary to avoid the warning and provide better estimation. 

Usually, high coverage would be a better choice to reliably detect much longer repeats (hundreds of repeats); otherwise, each estimation repeat count only has very few supporting reads, and RepeatHMM might filter the estimation at the final step, and just outputs something like [0, 0], indicating that there is not enough evidence to support the estimation.


## Other explanation of logging files

Some information below is for the DEBUG purpose, and they would not be in the output file with default settings.  

### INFO: ref_repeat 14 CTG +14

The information of the CTG repeat in the reference genome, and "+14" is the information provided by users or in pre-defined repeats.

### INFO: no_repeat_id_list 40

The number of reads whose repeats could not be detected by TRF (Tandem Repeat Finder).

### INFO: short_repeat_id_list 0

The number of reads whose repeats could be detected by TRF but the repeat regions are too short.

### INFO: repeat_id_list 14

The number of reads whose repeats could be detected by TRF.

### INFO: expRegionInLongRead= 10 4

The number of reads with longer repeat regions detected, followed by the number of reads without longer repeat regions or incorrectly aligned or partial alignment with the region of interest.

### INFO: nonRepeatinLongRead= 7 10

The number of reads with possibly short repeat regions detected, followed by the number of reads with incorrect alignment or partial alignment with the region of interest.

### INFO: tol_info=['', 0, 0]

The tolerant repeat patterns. For example, CTT could be considered to well match with CTG in the error correction process (not  in the HMM process), although the last alignment are T vs G. Here, not tolerant repeat patterns were specified. 

### INFO:  allocr:20:1, 22:1, 23:2, 24:7, 25:76, 26:625, 27:520, 28:51, 29:6, 30:3, 31:3, 32:1, 35:2, 36:1, 37:2, 39:2, 40:3, 41:2, 42:3, 43:5, 45:4, 46:6, 47:4, 48:5, 49:3, 50:6, 51:4, 52:6, 53:2, 54:6, 55:4, 56:1, 57:2, 58:4, 59:4, 60:4, 61:3, 62:2, 63:3, 64:2, 65:2, 66:5, 67:11, 68:24, 69:35, 70:43, 71:35, 72:27, 73:7, 74:4, 75:5, 78:1, 83:1

The dictionary of the estimation repeat count. In 'X:Y', 'X' is the estimation repeat count, and 'Y' is the number of supporting reads. 




