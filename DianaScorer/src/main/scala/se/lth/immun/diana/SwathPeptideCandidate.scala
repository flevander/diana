package se.lth.immun.diana

import se.lth.immun.mzml.ghost.XChromatogramGroup
import se.lth.immun.anubis.ReferencePrecursor


object SwathPeptideCandidate {
	
	import SwathPeakCandidate._
	
	def apply(
			rp:ReferencePrecursor, 
			cg:XChromatogramGroup, 
			c:Carrier,
			pScoreFunc:Carrier => Double,
			n:Int = 1
	) = {
		val x 		= new SwathPeptideCandidate(rp)
		
		x.fragmentRankAllRatioProb 		= c.fragmentRankAllRatioProb 
		x.fragmentRankPcsRatioProb 		= c.fragmentRankPcsRatioProb 
		x.fragmentMarkovAllRatioProb 	= c.fragmentMarkovAllRatioProb 
		x.fragmentMarkovPcsRatioProb 	= c.fragmentMarkovPcsRatioProb 
		x.fragmentCorrScore 			= c.fragmentCorrScore
		
		x.isotopeRankAllRatioProb 	= c.isotopeRankAllRatioProb 
		x.isotopeRankPcsRatioProb 	= c.isotopeRankPcsRatioProb 
		x.isotopeMarkovAllRatioProb = c.isotopeMarkovAllRatioProb 
		x.isotopeMarkovPcsRatioProb = c.isotopeMarkovPcsRatioProb 
		x.isotopeCorrScore			= c.isotopeCorrScore

		x.rtProb				= c.rtProb
		x.pscore 	= pScoreFunc(c)
		x.correctedAreas 	= c.fragmentEstimation.areas
		x.rawAreas 	= cg.chromatograms.map(_.intensities.slice(c.g.istart, c.g.iend).sum)
		x.missing 	= false
		
		val t 		= cg.chromatograms.head.times
		x.rtStart 	= t(c.g.istart)
		x.rtEnd 	= t(c.g.iend)
		
		var max 	= 0.0
		var imax 	= -1
		for (chrom <- cg.chromatograms)
			for (i <- c.g.istart until c.g.iend)
				if (chrom.intensities(i) >= max) {
					max = chrom.intensities(i)
					imax = i
				}
		x.rtApex 		= t(imax)
		x.maxIntensity 	= max
		x.maxEstimate 	= c.fragmentEstimation.estimateApex
		x.alternatives 	= n
		
		x
	}
}


class SwathPeptideCandidate(
		var reference:ReferencePrecursor
) {
	
	var missing 	= true
	var fragmentRankAllRatioProb 	= 1.0
	var fragmentRankPcsRatioProb 	= 1.0
	var fragmentMarkovAllRatioProb 	= 1.0
	var fragmentMarkovPcsRatioProb 	= 1.0
	var fragmentCorrScore 	= 1.0
	
	var isotopeRankAllRatioProb 	= 1.0
	var isotopeRankPcsRatioProb 	= 1.0
	var isotopeMarkovAllRatioProb 	= 1.0
	var isotopeMarkovPcsRatioProb 	= 1.0
	var isotopeCorrScore	= 0.0
	
	var rtProb				= 0.0
	var pscore 		= 1.0
	var qvalue 		= -1.0
	var correctedAreas 		= Seq[Double]()
	var rawAreas 		= Seq[Double]()
	var maxIntensity 	= -1.0
	var maxEstimate		= -1.0
	var alternatives	= 0 // total number of candidates for this assay
	
	var rtStart 	= -1.0
	var rtApex 		= -1.0
	var rtEnd 		= -1.0
	
	def correctedArea = if (correctedAreas == null || correctedAreas.isEmpty) 0.0 else correctedAreas.sum
	def rawArea = if (rawAreas == null || rawAreas.isEmpty) 0.0 else rawAreas.sum
}
