#!/usr/bin/python

import numpy
import scipy
import scipy.stats
import scipy.optimize
import matplotlib as mpl
mpl.use('Agg')
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
	dit = iter(d.iteritems())
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
	

def findIRtModel(df, outBase, minPeps):
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
                                figPath = "%s_round%d_irt_mapping.png" % (outBase, iter)
                                #print "figPath:", figPath
				fig.savefig(figPath)
		except Exception as e:
			print "error drawing figure:", e
		
		return (p, corr*corr)
	
	
	ok 	= df[df['q-value'] < 0.05]
	print " iter    n       corr2     slope    intersect"
	if len(ok) < minPeps:
		return {"intercept":0, "slope":0, "std":0, "n":len(ok), "corr2":0}
	
	p, c2 	= fitAndPlot(ok['rtApex'], ok['rtApexAssay'], 1)
	dt 	= ok['rtApexAssay'] - (p[0]*ok['rtApex'] + p[1])
	
	ok2 	= ok[numpy.abs(dt - dt.mean()) < 2*dt.std()]
	
	p, c2 	= fitAndPlot(ok2['rtApex'], ok2['rtApexAssay'], 2)
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
	
	return {"intercept":norm[0] + p[-1], "slope":p[-2], "std":norm[1], "n":len(ok), "corr2":c2}






usage = """usage:
> irtMap.py [--diana-dir=X] [--min-peps=X(def 5)] real.csv decoy.csv [outFile]"""

outFile = ""
minPeps = 5

def readArgs(args):
	if args[0].lower().startswith("--diana-dir="):
		readArgs(args[1:])
	elif args[0].lower().startswith("--min-peps="):
		global minPeps
		minPeps = int(args[0].split("=")[1])
		readArgs(args[1:])
	elif len(args) < 2:
		print "too few arguments!"
		print usage
		exit(1)
	elif len(args) > 3:
		print "too many arguments!"
		print usage
		exit(1)
	else:
		global realFile
		global decoyFile
		global outFile
		realFile 		= args[0]
		decoyFile 	= args[1]
		if len(args) == 3:
			outFile = args[2]

if len(sys.argv) < 3:
	print usage
	exit(1)
else:	
	readArgs(sys.argv[1:])


real = read_csv(realFile, sep="\t")
decoy = read_csv(decoyFile, sep="\t")

real['pscore'] = 1 - real.fragmentMarkovAllRatioProb
decoy['pscore'] = 1 - decoy.fragmentMarkovAllRatioProb

real = calculateQvalues(real, decoy)

irtModel = findIRtModel(real, withoutExt(outFile), minPeps)

if irtModel["n"] < minPeps:
	print "ERROR: couldn't detect enough rt-peptides to create mapping. Only found %d of %d." % (irtModel["n"], minPeps) 
	exit(1)

if outFile == "":
	print "rt peps w. fdr < 0.01", len(real[real['q-value'] < 0.01])
	print "slope: ", irtModel['slope']
	print "intercept: ", irtModel['intercept']
	print "std: ", irtModel['std']
	print "corr2: ", irtModel['corr2']
else:
	out = open(outFile, "w")
	out.write("rt peps w. fdr < 0.01: %d\n" % len(real[real['q-value'] < 0.01]))
	out.write("slope: %f\n" % irtModel['slope'])
	out.write("intercept: %f\n" % irtModel['intercept'])
	out.write("std: %f\n" % irtModel['std'])
	out.write("corr2: %f\n" % irtModel['corr2'])
	out.close()


