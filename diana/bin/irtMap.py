#!/usr/bin/python

import numpy
import scipy
import scipy.stats
import scipy.optimize
import matplotlib.pyplot as plt
from pandas import *
import matplotlib.mlab as mlab
import sklearn
import sys
from dianaLib import *


def calculateQvalues(real, decoy):
	real['q-value'] = 2.0
	rmax = real.groupby(['peptideSequence', 'q1'])['pscore'].transform(np.max)
	r = real[real.pscore == rmax].sort_index(by='pscore', ascending=False)
	d = decoy.groupby(['peptideSequence', 'q1'])['pscore'].aggregate(np.max).order(ascending=False)
	ireal 	= 0
	idecoy 	= 0
	rit = r.iterrows()
	dit = d.iteritems()
	ixr, xr = rit.next()
	ixd, xd = dit.next()
	gdratio = float(len(r)) / len(d)
	maxq 	= 1.0 / (len(r) * len(d))
	
	try:
		while True:
			if (xr['pscore'] > xd):
				maxq = max(
						maxq, 
						(float(idecoy) / (ireal+1)) * gdratio
						)
				x = r.set_value(ixr, 'q-value', maxq)
				ixr, xr = rit.next()
				ireal += 1
			else:
				ixd, xd = dit.next()
				idecoy += 1
	except StopIteration, e:
		pass #print e
	
	#r.describe()
	
	try:
		for ixr, xr in rit:
			x = r.set_value(ixr, 'q-value', 1.0)
	except:
		pass
	
	real['q-value'] = r['q-value']
	real['q-value'] = real['q-value'].fillna(2.0)
	
	return real
	

def findIRtModel(df, outBase):
	def fitAndPlot(x, y, iter):
		p = numpy.polyfit(x, y, 1)
		def fit(x):
			return p[1] + p[0]*x
		
		lrx = numpy.array([x.min(), x.max()])
		corr = x.corr(y)
		title = "corr2: %.2f  intersect: %.2f   slope: %.2f" % (corr*corr, p[1], p[0])
		
		print "% 3d     % 5d   %.7f %.7f %.7f" % (iter, len(x), corr*corr, p[0], p[1])
		
		try:
			fig = plt.figure()
			plt.title(title)
			plt.plot(lrx, fit(lrx), 'k-', x, y, 'bo', alpha=0.2, lw=3)#; plt.show()
			if outBase != "":
				fig.savefig("%s_round%d_irt_mapping.png" % (outBase, iter))
		except:
			pass
		
		return p
	
	
	ok 	= df[df['q-value'] < 0.05]
	print " iter    n       corr2     slope    intersect"
	if len(ok) < 10:
		return {"intercept":0, "slope":0, "std":0}
	
	p 	= fitAndPlot(ok['rtApex'], ok['rtApexAssay'], 1)
	dt 	= ok['rtApexAssay'] - (p[0]*ok['rtApex'] + p[1])
	
	ok2 	= ok[numpy.abs(dt - dt.mean()) < 2*dt.std()]
	
	p 	= fitAndPlot(ok2['rtApex'], ok2['rtApexAssay'], 2)
	dt 	= ok['rtApexAssay'] - (p[0]*ok['rtApex'] + p[1])
	
	dt.index = range(dt.size)
	n, bins = numpy.histogram(dt, 200, normed=True)
	# n, bins, patches = ax.hist(dt, 200, normed=1, facecolor='green', alpha=0.75)
	ns = Series(n)
	bincenters = 0.5 * (bins[1:] + bins[:-1])
	
	def error(x):     
		mu = x[0]    
		sigma = x[1]
		y = mlab.normpdf( bincenters, mu, sigma)
		ys = Series(y)
		return ((ns - ys)*(ns - ys)).sum()
	
	norm = scipy.optimize.fmin(error, [dt.mean(), dt.std()], disp=False)
	y = mlab.normpdf( bincenters, norm[0], norm[1])
	
	try:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.hist(df, 200, normed=1, facecolor='green', alpha=0.75)
		l = ax.plot(bincenters, y, 'ro', linewidth=1)
	
		if outBase != "":
			fig.savefig("%s_irt_mapping_error.png" % outBase)
		
	except:
		pass
	
	return {"intercept":norm[0] + p[-1], "slope":p[-2], "std":norm[1]}


if len(sys.argv) < 3:
	print "irtMap.py real.csv decoy.csv [outFile]"
	exit(1)

outFile = ""
if len(sys.argv) == 4:
	outFile = sys.argv[3]

realFile 	= sys.argv[1]
decoyFile 	= sys.argv[2]

real = read_csv(realFile, sep="\t")
decoy = read_csv(decoyFile, sep="\t")

real['pscore'] = 1 - real.fragmentMarkovAllRatioProb
decoy['pscore'] = 1 - decoy.fragmentMarkovAllRatioProb
real = calculateQvalues(real, decoy)

irtModel = findIRtModel(real, withoutExt(outFile))

if outFile == "":
	print "rt peps w. fdr < 0.01", len(real[real['q-value'] < 0.01])
	print "slope: ", irtModel['slope']
	print "std: ", irtModel['std']
	print "intercept: ", irtModel['intercept']
else:
	out = open(outFile, "w")
	out.write("rt peps w. fdr < 0.01: %d\n" % len(real[real['q-value'] < 0.01]))
	out.write("slope: %f\n" % irtModel['slope'])
	out.write("intercept: %f\n" % irtModel['intercept'])
	out.write("std: %f\n" % irtModel['std'])
	out.close()


