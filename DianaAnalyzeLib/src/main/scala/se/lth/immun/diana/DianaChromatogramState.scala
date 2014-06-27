package se.lth.immun.diana

import se.lth.immun.math.Ratios
import se.lth.immun.math.Matrix
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.distribution.NormalDistribution
import se.lth.immun.math.Stats
import DianaPeakCandidate._
import DianaInput._
import se.lth.immun.markov.MarkovChainDistribution




object DianaChromatogramState {
	class Stats {
		var ratios			:Ratios 								= null
		var variance												= -1.0
		var correlations	:Array[Array[Double]] 					= null
		var markovDists		:Array[Array[MarkovChainDistribution]] 	= null
		
		def clear = {
			ratios = null
			variance = -1.0
			correlations = null
			markovDists = null
		}
		
		def calculate(
				pcMaxWidth:Int, 
				numChromatograms:Int, 
				targetRatioTable:Array[Array[Ratio]]
		) = {
			// Markov dists
			markovDists = DianaUtil.calculateMarkovDists(
									numChromatograms, 
									ratios, 
									targetRatioTable, 
									pcMaxWidth)
			
									
			// Variance
			variance = 0.0
			Ratios.iterate(numChromatograms, (ri, ui, di) => {
				var r = ratios.getRatio(ri)
							.filter(d => d > 0.0 && d < Double.PositiveInfinity)
				variance += StatUtils.variance(r)
			})
			variance /= numChromatograms * (numChromatograms - 1) / 2
			
			// Correlations
			val rL			= ratios.length
			correlations 	= Matrix.get2d[Double](rL)
			Ratios.iterate(rL, (ri, ui, di) => {
				var rr = ratios.getRatio(ui).zip(ratios.getRatio(di))
							.filter(t => 
										t._1 > 0.0 && t._1 < Double.PositiveInfinity
									&&	t._2 > 0.0 && t._2 < Double.PositiveInfinity
								).unzip
				correlations(ui)(di) = Stats.pearsonCorrelation(rr._1.toArray, rr._2.toArray)
			})
		}
	}
}




class DianaChromatogramState(
		var _upperBound:Double,
		var rtIntersect:Double 	= Double.NaN, //expected empirical rt mean in minutes
		var rtSlope:Double 		= Double.NaN, //expected empirical rt slope in minutes
		var rtStd:Double 		= Double.NaN, //expected empirical rt std in minutes
		var binSize:Int 	= DianaUtil.DEFAULT_BIN_SIZE
) {
	import DianaChromatogramState._
	
	var _lowerBound = 1 / _upperBound
	def upperBound_=(d:Double) = { _upperBound = d; _lowerBound = 1 / d; }
	def lowerBound_=(d:Double) = { _lowerBound = d; _upperBound = 1 / d; }
	def upperBound = _upperBound
	def lowerBound = _lowerBound
	
	var numChromatograms	:Int 							= 0
	var targetRatioTable	:Array[Array[Ratio]] 			= null
	var expectedRt			:Double							= Double.NaN // in minutes
	var times				:Array[Double] 					= null
	var npcgs				:Int							= 0
	var ratioOkDists		:Array[Array[Array[Double]]]	= null
	var statsAll	= new Stats
	var statsPcs	= new Stats

	var rtNormalDistribution:NormalDistribution = null 
	
	
	def setChromatogram(
			numChromatograms:Int,
			ratios:Ratios,
			targetRatioTable:Array[Array[Ratio]],
			expectedRt:Double = Double.NaN,
			times:Array[Double] = null
	) = {
		this.numChromatograms 	= numChromatograms
		this.targetRatioTable 	= targetRatioTable
		this.expectedRt 		= expectedRt
		this.times 				= times
		statsAll.clear
		statsPcs.clear
		statsAll.ratios = ratios
		ratioOkDists	= null
		if (rtMapUsed)
			rtNormalDistribution = new NormalDistribution(rtSlope * expectedRt + rtIntersect, rtStd)
	}
	
	
	
	def isNaN = java.lang.Double.isNaN _
	def rtMapUsed = {
		!(isNaN(rtIntersect) || isNaN(rtSlope)  || isNaN(rtStd))
	}
	
	
	
	def rtProb(it:Int) = {
		if (!rtMapUsed || isNaN(expectedRt) || times == null)
			1.0
		else
			rtNormalDistribution.density(times(it))
	}
	
	
	
	def ratioOkProb(qouta:Double, ui:Int, di:Int) = {
		ratioOkDists(ui)(di).count(_ >= qouta).toDouble / npcgs
	}
	
	
	
	def calculateChromatogramStatistics(pcs:Seq[(PCGroup, GroupValidation)]) = {
		
		val maxWidth 	= pcs.map(_._1.iwidth).max
		npcgs			= pcs.length
		ratioOkDists	= Matrix.get2d[Array[Double]](numChromatograms)
		Ratios.iterate(numChromatograms, (ri, ui, di) => {
			ratioOkDists(ui)(di) = pcs.map(t => t._1.getQuota(ui, di, t._2)).sorted.toArray
		})
		
		
		val intervals = pcs.map(t => (t._1.istart, t._1.iend))
		val iL	= intervals.map(t => t._2 - t._1).sum
		statsPcs.ratios = new Ratios(numChromatograms)
		Ratios.iterate(numChromatograms, (ri, ui, di) => {
			statsPcs.ratios.addRatio(slices(statsAll.ratios.getRatio(ri), intervals, iL), ui, di)
		})
		
		
		statsAll.calculate(maxWidth, numChromatograms, targetRatioTable)
		statsPcs.calculate(maxWidth, numChromatograms, targetRatioTable)
	}
	
	
	
								
	def slices(x:Array[Double], intervals:Seq[(Int, Int)], iLength:Int = -1) = {
		var iL = iLength
		if (iL < 0) iL = intervals.map(t => t._2 - t._1).sum
		val ret = new Array[Double](iL)
		var ir = 0
		for (ival <- intervals) 
			for (ix <- ival._1 until ival._2) {
				ret(ir) = x(ix)
				ir += 1
			}
		ret
	}
}
