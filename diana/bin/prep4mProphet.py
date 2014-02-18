#!/usr/bin/python

import subprocess
import sys
from pandas import *
import scipy.stats




def calcIRtProb(df, irtKs):
	irtProb = lambda x:(1-scipy.stats.norm.cdf(abs(x)))*2
	df['measuredIRT'] = irtKs["intercept"] + df['rtApex']*irtKs["slope"]
	df['irtDev'] = df['rtApexAssay'] - df['measuredIRT']
	irtStat = df['irtDev'] / irtKs["std"]
	df['irtProb'] = irtStat.apply(irtProb)




def setupDF(df, isDecoy, pCutoff=0.99):
	df 		= df.fillna(1.0)
	fcs 	= df.fragmentCorrScore
	fcs[fcs > 1.0] = 1.0
	fcs[fcs < -1.0] = -1.0
	df['inArea']		= np.log(df.area + 1.0) * 10.0 # because mProphet LDA complains
	df['var_inRt'] 		= df.irtProb * 100 # because mProphet LDA complains
	df['var_irtAbsDev'] = -df.irtDev.abs()
	df['var_logRt'] 	= np.log(df.irtProb + 0.0000000001)
	df['main_var_logFMPRP'] 	= np.log(df.fragmentMarkovPcsRatioProb + 0.0000000001)
	df['main_var_inFMPRP'] 		= 1.0 - df.fragmentMarkovPcsRatioProb
	df['var_inFCS'] 	= df.fragmentCorrScore
	df['var_logIMPRP'] 	= np.log(df.isotopeMarkovPcsRatioProb + 0.0000000001)
	df['var_inIMPRP'] 	= 1.0 - df.isotopeMarkovPcsRatioProb
	df['var_inICS'] 	= df.isotopeCorrScore
	df['var_inAbsRt'] 	= df.rtApex / 200
	df['decoy']			= isDecoy
	df['transition_group_record'] = df['protein'] + '|' + df['peptideSequence'] + '|' + df['charge'].apply(str)
#	df['logArea'] = df['logArea'] - df['logArea'].mean()) / df['logArea'].std()
	return df




def setup(irtFile):
	real 	= read_csv(realFile, sep="\t")
	decoy 	= read_csv(decoyFile, sep="\t")
	
	irtKs = {}
	irt = open(irtFile, 'r')
	for line in irt:
		if ':' in line:
			t = line.split(":")
			irtKs[t[0].strip()] = float(t[1].strip())
	
	print "slope: ", 	float(irtKs["slope"])
	print "icept: ", 	float(irtKs["intercept"])
	print "std: ", 		float(irtKs["std"] )
	
	calcIRtProb(real, irtKs)
	calcIRtProb(decoy, irtKs)
	
	return setupDF(real, False), setupDF(decoy, True)





if len(sys.argv) < 4:
	print "extract_eval.py real.csv decoy.csv irt.csv"
	exit(1)

realFile 	= sys.argv[-3]
decoyFile 	= sys.argv[-2]
irtCsvFile	= sys.argv[-1]
irtFile 	= "irtmap"

subprocess.call(["/mnt/africa/code/scripts/irtMap.py", irtCsvFile, decoyFile, irtFile])

real, decoy = setup(irtFile)

merged = concat([real, decoy], ignore_index=True)

merged.to_csv(
	"%s.to.mProph" % realFile, 
	index=False, 
	sep='\t',
	cols=['transition_group_record', 'decoy', 'main_var_inFMPRP', 
			'var_irtAbsDev', 'var_inFCS', 'var_inIMPRP', 'var_inICS']
	)
																	
