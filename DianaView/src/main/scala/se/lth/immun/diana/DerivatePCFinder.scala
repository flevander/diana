package se.lth.immun.diana

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.distribution.BinomialDistribution
import org.apache.commons.math3.distribution.ChiSquaredDistribution

import se.lth.immun.math.Matrix
import se.lth.immun.math.Ratios
import se.lth.immun.anubis.Ratio
import se.lth.immun.markov.MarkovChainDistribution


object DerivatePCFinder {

	var UPPER_BOUND = 1.5
    var LOWER_BOUND = 1.0 / UPPER_BOUND
    
    
	class PC(
			val istart:Int,
			val iapex:Int,
			val iend:Int
	) {
		var icurve = -1
	}
	
	
	
	class PCGroup {
		var pcs 	= new ArrayBuffer[PC]
		var pvalue 	= -1.0
		var pvalueTable:Array[Array[Double]] 	= null
		var validations:Array[Validation]		= null
		var ratioValidations:Array[Array[Validation]] = null
		
		override def toString = istart + "-"+ iend
			
		def iend 	= StatUtils.percentile(pcs.map(_.iend.toDouble).toArray, 50).toInt //pcs.map(_.iend).max
		def istart 	= StatUtils.percentile(pcs.map(_.istart.toDouble).toArray, 50).toInt //pcs.map(_.istart).min
		def iwidth 	= iend - istart
		
		var nistart = 0
		var niend = 0
		def niwidth = niend - nistart
	}
	
	
	
	class Validation(
		var icurve:Int,
		width:Int
	) {
		var ok 		= new Array[Boolean](width)
		var flag 	= false
		for (i <- 0 until width) ok(i) = false
	}
	
	
	
	
	
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
			targetRatioTable:Array[Array[Ratio]]
	):Array[Array[Double]] = {
		var table 	= Matrix.get2d[Double](numChromatograms)
		Ratios.iterate(numChromatograms, (ri, ui, di) => {
			var r = ratios.getRatio(ui, di)
			var target = targetRatioTable(ui)(di).mean
			//markovPValues(r, target)
			var c = r.count(d => 	d > 0.0
								&&	d < Double.PositiveInfinity
								&&	d < target * UPPER_BOUND
	            				&&	d > target * LOWER_BOUND)
	        table(ui)(di) = math.max(c, 1.0) / r.length
		})
		return table
	}
	
	
	
	def markovPValues(r:Array[Double], target:Double) = {
		var counts = Array(Array(0, 0), Array(0, 0))
		for (i <- 1 until r.length) {
			val r0 = r(i-1)
			val r1 = r(i)
			var t0 = if (	r0 > 0.0
					&&	r0 < Double.PositiveInfinity
					&&	r0 < target * UPPER_BOUND
					&&	r0 > target * LOWER_BOUND) 1 else 0
			var t1 = if (	r1 > 0.0
					&&	r1 < Double.PositiveInfinity
					&&	r1 < target * UPPER_BOUND
					&&	r1 > target * LOWER_BOUND) 1 else 0
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
			n:Int
	):MarkovChainDistribution = {
		var counts = Array(Array(0, 0), Array(0, 0))
		for (i <- 1 until r.length) {
			val r0 = r(i-1)
			val r1 = r(i)
			var t0 = if (	r0 > 0.0
					&&	r0 < Double.PositiveInfinity
					&&	r0 < target * UPPER_BOUND
					&&	r0 > target * LOWER_BOUND) 1 else 0
			var t1 = if (	r1 > 0.0
					&&	r1 < Double.PositiveInfinity
					&&	r1 < target * UPPER_BOUND
					&&	r1 > target * LOWER_BOUND) 1 else 0
			counts(t0)(t1) += 1
		}
		val sums = counts.map(_.sum)
		val p1 = (counts(0)(1)+counts(1)(1)) / sums.sum.toDouble
		val p11 = if (sums(1) == 0) 0.0 else counts(1)(1) / sums(1).toDouble
		val p22 = if (sums(0) == 0) 0.0 else counts(0)(0) / sums(0).toDouble
		return new MarkovChainDistribution(n, p1, p11, p22)
	}
}



class DerivatePCFinder(
		var _upperBound:Double,
		var binSize:Int
) {
	var _lowerBound = 1 / _upperBound
	def upperBound_=(d:Double) = { _upperBound = d; _lowerBound = 1 / d; }
	def lowerBound_=(d:Double) = { _lowerBound = d; _upperBound = 1 / d; }
	def upperBound = _upperBound
	def lowerBound = _lowerBound
	
	
	import DerivatePCFinder._
	
	var numChromatograms:Int 	= 0
	var ratios:Ratios 			= null
	var targetRatioTable:Array[Array[Ratio]] = null
	var markovDists:Array[Array[MarkovChainDistribution]] = null
	var variance			= -1.0
	var correlations:Array[Array[Double]] = null
	
	
	
	def setChromatogram(
			numChromatograms:Int,
			ratios:Ratios,
			targetRatioTable:Array[Array[Ratio]]
	) = {
		this.numChromatograms 	= numChromatograms
		this.ratios 			= ratios
		this.targetRatioTable 	= targetRatioTable
		markovDists 	= null
		variance		= -1.0
		correlations 	= null
	}
	
	
	
	def calculateChromatogramStatistics(maxPCWidth:Int) = {
		markovDists	= calculateMarkovDists(numChromatograms, 
								ratios, targetRatioTable, maxPCWidth)
		
								
		val rL			= ratios.length
		val pearson 	= new PearsonsCorrelation
		correlations 	= Matrix.get2d[Double](rL)
		variance 		= 0.0
		
		Ratios.iterate(numChromatograms, (ri, ui, di) => {
			var r = ratios.getRatio(ri)
						.filter(d => d > 0.0 && d < Double.PositiveInfinity)
			variance += StatUtils.variance(r)
		})
		variance /= 45
		
		Ratios.iterate(rL, (ri, ui, di) => {
			var rr = ratios.getRatio(ui).zip(ratios.getRatio(di))
						.filter(t => 
									t._1 > 0.0 && t._1 < Double.PositiveInfinity
								&&	t._2 > 0.0 && t._2 < Double.PositiveInfinity
							).unzip
			correlations(ui)(di) = pearson.correlation(rr._1.toArray, rr._2.toArray)
		})
	}
	
	
	def findPCs(
			y:Seq[Double],
			dy:Seq[Double],
			ddy:Seq[Double],
			baseline:Seq[Double],
			_binSize:Int = 0
	):Seq[PC] = {
		var ret = new ArrayBuffer[PC]
		
		var istart 	= -1
		var iapex 	= -1
		val bs		= if (_binSize > 0) _binSize else binSize
		for (i <- 0 until dy.length-1) {
			if (math.signum(dy(i)) != math.signum(dy(i+1))) {
				if (ddy(i) > 0) { // local min
					if (iapex >= 0 && istart >= 0) {
						var max 	= y(iapex)
						var starty 	= y(istart)
						var endy	= y(i+1)
						var bl = baseline(math.max(0, math.min(baseline.length-1, iapex - bs / 2)))
						if (max > 2*math.max(starty, endy) && max > 2 * bl)
							ret += new PC(istart, iapex, i+1)
					}
						
					istart = i+1
					iapex = -1
				} else { // local max
					iapex = i+1
				}
			}
		}
		
		return ret
	}
	
	
	
	def groupPCs(
			pcs:Seq[Seq[PC]]
	):Seq[PCGroup] = {
		for (i <- 0 until pcs.length)
			pcs(i).foreach(_.icurve = i)
			
		var sorted 	= pcs.flatten.sortBy(_.iapex)
		var edges	= findEdges(sorted)
		
		var L				= sorted.length
        var groupRefList 	= new Array[PCGroup](L)
        var groupList 		= new ArrayBuffer[PCGroup]

        for (i <- 0 until L) {
            for (e <- edges(i))
                if (groupRefList(e) != null) groupRefList(i) = groupRefList(e)
            
            if (groupRefList(i) == null) {
                var pc = new PCGroup
                groupList += pc
                groupRefList(i) = pc
            }
            for (e <- edges(i))
                if (groupRefList(e) == null) groupRefList(e) = groupRefList(i)
        }

        for (i <- 0 until L)
            groupRefList(i).pcs += sorted(i)
        
        return groupList
	}
	
	
	
	def findEdges(
			pcs:Seq[PC]
	):Array[List[Int]] = {
		var L 		= pcs.length
		var edges 	= new Array[List[Int]](L)
		
		for (i <- 0 until L) {
            edges(i) = Nil
            var k = i + 1
            var li = pcs(i)
            var ki = li
            while ({
            	if (k < L) {
            		ki = pcs(k)
            		li.iapex > ki.istart && ki.iapex < li.iend
            	} else false
            }) {
            	var ui = math.min(li.icurve, ki.icurve)
            	var di = math.max(li.icurve, ki.icurve)
            	var r 		= ratios.getRatio(ui, di)
            	var target 	= targetRatioTable(ui)(di).mean
            	var count 	= 0
            	for (ir <- ki.istart until li.iend+1) {
            		if (		r(ir) < target * UPPER_BOUND
            				&&	r(ir) > target * LOWER_BOUND )
            			count += 1
            	}
                if (2 * count >= li.iend - ki.istart)
                	edges(i) = edges(i) :+ k
                
                k += 1
            }
        }
		return edges
	}
	
	
	
	def validateGroup(
			group:PCGroup,
			smooths:Seq[(Array[Double], Seq[Double])]
	):PCGroup = {
		var istart 	= group.istart
		var iend 	= group.iend 
		var iwidth	= iend - istart+1
		var L		= targetRatioTable.length
		
		//var counts = Matrix.get2d[Int](L).map(_.map(i => -1))
		var handled = targetRatioTable.map(r => false)
		var valids = new Array[Validation](L)
		var rValids = Matrix.get2d[Validation](L).map(_.map(v => new Validation(0, iwidth)))
		for (i <- 0 until L)
			valids(i) = new Validation(i, iwidth)

		
		var q = new Queue[Validation]
		for (pc <- group.pcs) {
			q += valids(pc.icurve)
			for (ir <- math.max(istart, pc.istart) until math.min(iend, pc.iend)+1)
				valids(pc.icurve).ok(ir - istart) = true
		}
		
		while (!q.isEmpty) {
			var v = q.dequeue
			for (i <- 0 until L) {
				if (i == v.icurve)
					rValids(i)(i).flag = true // = v.ok.count(b => b)
				else if (!rValids(v.icurve)(i).flag) {
					var ui = math.min(i, v.icurve)
	            	var di = math.max(i, v.icurve)
	            	var r 		= ratios.getRatio(ui, di)
	            	var target 	= targetRatioTable(ui)(di).mean
	            	//var count 	= 0
	            	var y		= smooths(i)._1
	            	var bl		= smooths(i)._2
	            	var bL		= bl.length
	            	for (ir <- istart until iend+1) {
	            		if (		r(ir) < target * UPPER_BOUND
	            				&&	r(ir) > target * LOWER_BOUND 
	            				&&	v.ok(ir - istart)
	            				&& 	y(ir) > bl(math.max(0, math.min(bL-1, ir - binSize / 2)))) {
	            			//count += 1
	            			valids(i).ok(ir - istart) = true
	            			rValids(ui)(di).ok(ir - istart) = true
	            		} else
	            			rValids(ui)(di).ok(ir - istart) = false
	            	}
					rValids(v.icurve)(i).flag = true//counts(v.icurve)(i) = count
					rValids(i)(v.icurve).flag = true//counts(i)(v.icurve) = count
					if (rValids(i)(i).flag)
						q += valids(i)
				}
			}
		}
		
		var rvs 	= new ArrayBuffer[Array[Boolean]]
		//var minMax 	= new ArrayBuffer[(Int, Int)]
		Ratios.iterate(L, (ri, ui, di) => {
			var v = rValids(ui)(di)
			rvs += v.ok
			//if (v.ok.exists(b => b))
			//	minMax += new Tuple2(v.ok.indexWhere(b => b), v.ok.lastIndexWhere(b => b))
		})
		var rvsSummed 	= rvs.transpose.map(_.count(b => b))
		var vTotal 		= rvsSummed.sum
		var temp 		= 2.0 / rvs.length
		var ratio 		= vTotal.toDouble / (rvs.length * iwidth)
		
		var nistart 	= 0
		var niend 		= iwidth
		while (rvsSummed(nistart) * temp < ratio) nistart += 1
		while (rvsSummed(niend-1) * temp < ratio) niend -= 1
		
		group.nistart 	= nistart
		group.niend 	= niend
		group.validations 		= valids
		group.ratioValidations 	= rValids
		
		
		/*
		var nistart = JTMath.median(minMax.map(_._1.toDouble).toArray).toInt
		var niend = JTMath.median(minMax.map(_._2.toDouble).toArray).toInt + 1
		var niwidth = niend - nistart
		mrm.data.Ratios.iterate(L, (ri, ui, di) => {
			ratioValidations(ui)(di).ok = ratioValidations(ui)(di).ok.slice(nistart, niend)
		})
		*/
		return group
	}
	
	
	
	def pvalueForValidatedGroup(
			group:PCGroup
	):PCGroup = {
		var pvalues 	= new Array[Double](ratios.length)
		var corr 		= new Array[Double](ratios.length)
		var L			= targetRatioTable.length
		group.pvalueTable = Matrix.get2d[Double](L)
		
		Ratios.iterate(L, (ri, ui, di) => {
			val ok = group.ratioValidations(ui)(di).ok
			val count = 
				math.max(
						0, 
						ok.slice(group.nistart, group.niend).count(b => b)-1
					)
			val pdf = markovDists(ui)(di).pdf(group.niwidth)
			pvalues(ri) = pdf.take(pdf.length - count).sum
			/*
			val dist = new BinomialDistribution(niwidth, ratioPValueTable(ui)(di))
			pvalues(ri) = 1-dist.cumulativeProbability(count)//counts(ui)(di)-1))
					*/
			if (pvalues(ri) == 0.0) {
				pvalues(ri) = java.lang.Double.MIN_NORMAL
				println("########## counts: "+count+
						"    "+markovDists(ui)(di).toString(group.niwidth))
			}
			group.pvalueTable(ui)(di) 	= pvalues(ri)
			corr(ri)					= correlations(ui)(di)
		})
		
		val E = 2.0*ratios.length
		val v = 4.0*ratios.length + 2.0 * corr.map(c => estimateCov(c, variance)).sum
		val f = (2.0 * E * E) / v
		val c = v / (2.0 * E)
		
		/*
		println("pvalues: "+pvalues.mkString(" "))
		println("corr: "+corr.mkString(" "))
		println("covs: "+corr.map(c => estimateCov(c, variance)).mkString(" "))
		println("E: "+E)
		println("v: "+v)
		println("f: "+f)
		println("c: "+c)
		println("fishers: "+fishers(pvalues))
		*/
		val chi = new ChiSquaredDistribution(f)
		group.pvalue = 1-chi.cumulativeProbability(fishers(pvalues) / c)
		return group
	}
}
