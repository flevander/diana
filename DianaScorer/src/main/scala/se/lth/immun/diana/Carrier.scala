package se.lth.immun.diana

import se.lth.immun.anubis.ReferencePrecursor
import DianaPeakCandidate._

class Carrier(
		var assay:DianaAssay,
		var g:PCGroup,
		nFragments:Int,
		nIsotopes:Int
) {
	var fragmentValids:GroupValidation = null
	var fragmentEstimation:GroupEstimation = null
	var fragmentRankAllRatioProb 	= 1.0
	var fragmentRankPcsRatioProb 	= 1.0
	var fragmentMarkovAllRatioProb 	= 1.0
	var fragmentMarkovPcsRatioProb 	= 1.0
	var fragmentCorrScore		= 0.0
	var fragmentPcs = g.pcs.filter(_.icurve < nFragments)
	
	
	var isotopeValids:GroupValidation = null
	var isotopeEstimation:GroupEstimation = null
	var isotopeRankAllRatioProb 	= 1.0
	var isotopeRankPcsRatioProb 	= 1.0
	var isotopeMarkovAllRatioProb 	= 1.0
	var isotopeMarkovPcsRatioProb 	= 1.0
	var isotopeCorrScore	= 0.0
	var isotopePcs = (0 until nIsotopes).map(i => {
			val pc = new PC(g.istart, (g.iend + g.istart)/2, g.iend)
			pc.icurve = i
			pc
		})
	/*var isotopePcs = g.pcs.filter(_.icurve >= nFragments).map(pc => {
			pc.icurve -= nFragments
			pc
		})
	*/
	var rtProb 		= 0.0
	var pEstimates	= 1.0
	
	
	def isRemotelyUnlikely(cutoff:Double) = 
		fragmentRankAllRatioProb < cutoff ||
		fragmentRankPcsRatioProb < cutoff ||
		fragmentMarkovAllRatioProb < cutoff ||
		fragmentMarkovPcsRatioProb < cutoff
}
