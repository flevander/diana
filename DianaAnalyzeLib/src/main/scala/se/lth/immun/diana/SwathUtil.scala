package se.lth.immun.diana

import se.lth.immun.math.Ratios
import se.lth.immun.anubis.Ratio

import se.lth.immun.math.Matrix
import se.lth.immun.markov.MarkovChainDistribution

object SwathUtil {

	val DEFAULT_UPPER_BOUND = 1.5
    val DEFAULT_LOWER_BOUND = 1.0 / DEFAULT_UPPER_BOUND
	val DEFAULT_BIN_SIZE = 20
    
    
	/**
	 * Approximates the covariance between two stochastic variables
	 * Xi and Xj using their correlation pij and their variance v
	 */
	def estimateCov(pij:Double, v:Double):Double = 
		return pij * (3.263 + pij * (0.710 + pij*0.027)) + 
				(0.727 + pij*(0.327 - pij*(0.768 - pij*0.331)))/v
	
	
	/**
	 * Calculates the fishers statistic of a number of p-values:
	 * F = -2 * sum( ln(p) )
	 */
	def fishers(pvalues:Seq[Double]):Double = 
		return -2 * pvalues.map(p => math.log(p)).sum
	
	
	
	
	def getPValueTable(
			numChromatograms:Int,
			ratios:Ratios,
			targetRatioTable:Array[Array[Ratio]],
			upperBound:Double = -1.0
	):Array[Array[Double]] = {
		val ub = if (upperBound >= 0) upperBound else DEFAULT_UPPER_BOUND
		val lb = if (upperBound >= 0) 1 / upperBound else DEFAULT_LOWER_BOUND
		var table 	= Matrix.get2d[Double](numChromatograms)
		Ratios.iterate(numChromatograms, (ri, ui, di) => {
			var r = ratios.getRatio(ui, di)
			var target = targetRatioTable(ui)(di).mean
			//markovPValues(r, target)
			var c = r.count(d => 	d > 0.0
								&&	d < Double.PositiveInfinity
								&&	d < target * ub
	            				&&	d > target * lb)
	        table(ui)(di) = math.max(c, 1.0) / r.length
		})
		return table
	}
	
	
	
	def markovPValues(r:Array[Double], target:Double, upperBound:Double = -1.0) = {
		val ub = if (upperBound >= 0) upperBound else DEFAULT_UPPER_BOUND
		val lb = if (upperBound >= 0) 1 / upperBound else DEFAULT_LOWER_BOUND
		var counts = Array(Array(0, 0), Array(0, 0))
		for (i <- 1 until r.length) {
			val r0 = r(i-1)
			val r1 = r(i)
			var t0 = if (	r0 > 0.0
					&&	r0 < Double.PositiveInfinity
					&&	r0 < target * ub
					&&	r0 > target * lb) 1 else 0
			var t1 = if (	r1 > 0.0
					&&	r1 < Double.PositiveInfinity
					&&	r1 < target * ub
					&&	r1 > target * lb) 1 else 0
			counts(t0)(t1) += 1
		}
		println("    \t   TO")
		println("FROM\tok\tnot")
		println("not\t%d\t%d".format(counts(0)(1), counts(0)(0)))
		println("ok\t%d\t%d".format(counts(1)(1), counts(1)(0)))
	}
	
	
	
	def calculateMarkovDists(
			numChromatograms:Int,
			ratios:Ratios,
			targetRatioTable:Array[Array[Ratio]],
			n:Int
	):Array[Array[MarkovChainDistribution]] = {
		var table 	= Matrix.get2d[MarkovChainDistribution](numChromatograms)
		Ratios.iterate(numChromatograms, (ri, ui, di) => {
			var r 		= ratios.getRatio(ui, di)
			var target 	= targetRatioTable(ui)(di).mean
	        table(ui)(di) = calculateMarkovDist(r, target, n)
		})
		return table
	}
	
	
	
	def calculateMarkovDist(
			r:Array[Double], 
			target:Double, 
			n:Int,
			upperBound:Double = -1.0
	):MarkovChainDistribution = {
		val ub = if (upperBound >= 0) upperBound else DEFAULT_UPPER_BOUND
		val lb = if (upperBound >= 0) 1 / upperBound else DEFAULT_LOWER_BOUND
		var counts = Array(Array(0, 0), Array(0, 0))
		for (i <- 1 until r.length) {
			val r0 = r(i-1)
			val r1 = r(i)
			var t0 = if (	r0 > 0.0
					&&	r0 < Double.PositiveInfinity
					&&	r0 < target * ub
					&&	r0 > target * lb) 1 else 0
			var t1 = if (	r1 > 0.0
					&&	r1 < Double.PositiveInfinity
					&&	r1 < target * ub
					&&	r1 > target * lb) 1 else 0
			counts(t0)(t1) += 1
		}
		val sums = counts.map(_.sum)
		val p1 = (counts(0)(1)+counts(1)(1)) / sums.sum.toDouble
		val p11 = if (sums(1) == 0) 0.0 else counts(1)(1) / sums(1).toDouble
		val p22 = if (sums(0) == 0) 0.0 else counts(0)(0) / sums(0).toDouble
		return new MarkovChainDistribution(n, p1, p11, p22)
	}
}
