
import os;
import sys;

#from .myheader import *
from . import myheader

myDefaultsetting = {\
							"hg": 'hg38', \
							"hgfile": None, \
							"GapCorrection": 1, \
							"FlankLength": 30, \
							#"MatchInfo": None, \
							"outlog": myheader.M_WARNING, \
							"Tolerate": None, \
							"MinSup": 5, \
							"MaxRep": 10000, \
							"CompRep": '0', \
							"repeatName": None, \
							#"UserDefinedUniqID": None, \
							#"Patternfile": None, \
							"UserDefinedRepeat": myheader.UserDefinedRepeatdefault, \
							"SplitAndReAlign": 0, \
							"TRFOptions": "2_7_4_80_10_100", \
							"minTailSize": 70, \
							"minRepBWTSize": 70, \
							"RepeatTime": 5, \
							"BWAMEMOptions": "k8_W8_r7", \
							"hmm_insert_rate": 0.12, \
							"hmm_del_rate": 0.02, \
							"hmm_sub_rate": 0.02, \
							"SeqTech": None, \
							"transitionm": None, \
							"emissionm": None, \
							"region": None, \
							#"thread": 1, \
							"outFolder": 'align/', \
							"conted": 0, \
							"Onebamfile": None, \
							"SepbamfileTemp": None, \
						}

def getString(margs): 
	retstr = []
	defkeys = myDefaultsetting.keys(); 
	for df in defkeys:
		curv = eval('margs.'+df)
		if myDefaultsetting[df] == curv or curv==None:
			pass;
			#retstr.append('--'+df);
			#retstr.append(str(curv))
		else:
			retstr.append('--'+df)
			if isinstance(curv, str):
				#retstr.append('"'+str(curv)+'"')
				retstr.append(str(curv))
			else:
				retstr.append(str(curv))
	#print retstr
	return ' '.join(retstr)
