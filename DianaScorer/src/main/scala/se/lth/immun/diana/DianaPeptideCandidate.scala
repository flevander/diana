package se.lth.immun.diana

import se.lth.immun.mzml.ghost.XChromatogramGroup
import se.lth.immun.anubis.ReferencePrecursor


object DianaPeptideCandidate {
	
	import DianaPeakCandidate._
	
	def apply(
			assay:DianaAssay, 
			cg:XChromatogramGroup, 
			icg:Option[XChromatogramGroup], 
			c:Carrier,
			n:Int = 1
	) = {
		val x 		= new DianaPeptideCandidate(assay)
		
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
		x.correctedAreas 	= c.fragmentEstimation.areas
		x.rawAreas 	= cg.chromatograms.map(_.intensities.slice(c.g.istart, c.g.iend).sum)
		x.missing 	= false
		
		x.isotopeAreas = 
			icg match {
				case Some(x) => x.chromatograms.map(_.intensities.slice(c.g.istart, c.g.iend).sum)
				case None => Nil
			}
		
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


class DianaPeptideCandidate(
		var assay:DianaAssay
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
	var qvalue 		= -1.0
	var correctedAreas 	= Seq[Double]()
	var rawAreas 		= Seq[Double]()
	var isotopeAreas 	= Seq[Double]()
	var maxIntensity 	= -1.0
	var maxEstimate		= -1.0
	var alternatives	= 0 // total number of candidates for this assay
	
	var rtStart 	= -1.0
	var rtApex 		= -1.0
	var rtEnd 		= -1.0
	
	def safeSum(x:Seq[Double]) = if (x == null || x.isEmpty) 0.0 else x.sum
	def correctedArea 	= safeSum(correctedAreas)
	def rawArea 		= safeSum(rawAreas)
	def isotopeArea 	= safeSum(isotopeAreas)
}
