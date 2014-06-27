package se.lth.immun.diana

import DianaPeakCandidate._
import DianaSignalProcessor._
import DianaInput._
import se.lth.immun.math.Ratios

import se.lth.immun.math.Matrix
import se.lth.immun.math.Stats
import se.lth.immun.signal.Filter
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue

import org.apache.commons.math3.distribution.ChiSquaredDistribution
import org.apache.commons.math3.stat.StatUtils

class DianaPCEvaluator(
		var state:DianaChromatogramState,
		var minRatioValidity:Int = 2
) {

	
	
	
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
            		math.abs(li.iapex - ki.iapex) < 2
            	} else false
            }) {
                edges(i) = edges(i) :+ k
                k += 1
            }
        }
		return edges
	}
	
	
	
	
	
	
	def validateGroup(
			group:PCGroup,
			pcs:Seq[PC],
			smooths:Seq[SmoothAndBase]
	):GroupValidation = {
		var istart 	= group.istart
		var iend 	= group.iend 
		var iwidth	= group.iwidth
		var L		= state.targetRatioTable.length
		val ub 		= state.upperBound
		val lb 		= state.lowerBound
		
		var handled 	= state.targetRatioTable.map(r => false)
		var gv 			= GroupValidation(L, iwidth)
		
		var q = new Queue[Validation]
		for (pc <- pcs) {
			q += gv.valids(pc.icurve)
			for (ir <- math.max(istart, pc.istart) until math.min(iend, pc.iend))
				gv.valids(pc.icurve).ok(ir - istart) = true
		}
		
		while (!q.isEmpty) {
			var v = q.dequeue
			for (i <- 0 until L) {
				if (i == v.icurve)
					gv.rValids(i)(i).flag = true
				else if (!gv.rValids(v.icurve)(i).flag) {
					var ui = math.min(i, v.icurve)
	            	var di = math.max(i, v.icurve)
	            	var r 		= state.statsAll.ratios.getRatio(ui, di)
	            	var target 	= state.targetRatioTable(ui)(di).ratio
	            	var y		= smooths(i).smooth
	            	var bl		= smooths(i).base
	            	var bL		= bl.length
	            	for (ir <- istart until iend) {
	            		if (		r(ir) < target * ub
	            				&&	r(ir) > target * lb
	            				&&	v.ok(ir - istart)
	            				&& 	y(ir) > bl(math.max(0, math.min(bL-1, ir - state.binSize / 2)))) {
	            			gv.valids(i).ok(ir - istart) = true
	            			gv.rValids(ui)(di).ok(ir - istart) = true
	            		} else
	            			gv.rValids(ui)(di).ok(ir - istart) = false
	            	}
					gv.rValids(v.icurve)(i).flag = true
					gv.rValids(i)(v.icurve).flag = true
					if (gv.rValids(i)(i).flag)
						q += gv.valids(i)
				}
			}
		}
		
		for (i <- 0 until L) {
			var counts = new Array[Int](iwidth)
			for (j <- 0 until L) {
				if (i != j) {
					var ui = math.min(i, j)
		            var di = math.max(i, j)
		            var ok = gv.rValids(ui)(di).ok
					for (k <- 0 until iwidth) {
			            if (ok(k))
			            	counts(k) += 1
					}
				}
			}
			gv.valids(i).ok = counts.map(_ >= minRatioValidity)
		}
		
		var rvs 	= new ArrayBuffer[Array[Boolean]]
		Ratios.iterate(L, (ri, ui, di) => {
			var v = gv.rValids(ui)(di)
			rvs += v.ok
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
		
		return gv
	}
	
	
	
	
	
	
	
	def nullRatioProb(
			group:PCGroup,
			gv:GroupValidation,
			stats:DianaChromatogramState.Stats,
			mode:String
	):Double = {
		
		var pvalues 	= new Array[Double](stats.ratios.length)
		//var corr 		= new Array[Double](stats.correlations.length)
		var covSum		= 0.0
		var L			= state.numChromatograms
		group.pvalueTable = Matrix.get2d[Double](L)
		
		if (mode == "rank")
			Ratios.iterate(L, (ri, ui, di) => {
				pvalues(ri) = state.ratioOkProb(group.getQuota(ui, di, gv), ui, di)
	
				if (pvalues(ri) == 0.0) {
					pvalues(ri) = java.lang.Double.MIN_NORMAL
					println("########## qouta: "+group.getQuota(ui, di, gv))
				}
				group.pvalueTable(ui)(di) 	= pvalues(ri)
				//corr(ri)					= stats.correlations(ui)(di)
			})
		else if (mode == "markov")
			Ratios.iterate(L, (ri, ui, di) => {
				val ok = gv.rValids(ui)(di).ok
				val count = 
					math.max(
							0, 
							ok.slice(group.nistart, group.niend).count(b => b)-1
						)
				val pdf = stats.markovDists(ui)(di).pdf(group.niwidth)
				pvalues(ri) = pdf.take(pdf.length - count).sum
	
				if (pvalues(ri) == 0.0) {
					pvalues(ri) = java.lang.Double.MIN_NORMAL
					println("########## counts: "+count+
							"    "+stats.markovDists(ui)(di).toString(group.niwidth))
				}
				group.pvalueTable(ui)(di) 	= pvalues(ri)
				//corr(ri)					= stats.correlations(ui)(di)
			})
		else
			throw new Exception("Unknown DianaPCEvaluation.nullRatioProb() mode '"+mode+"'")
		
		Ratios.iterate(stats.correlations.length, (ri, ui, di) => 
				covSum += DianaUtil.estimateCov(stats.correlations(ui)(di), stats.variance)
			)
		
		val E = 2.0 * stats.ratios.length
		val v = 4.0 * stats.ratios.length + 2.0 * covSum
		val f = (2.0 * E * E) / v
		val c = v / (2.0 * E)
		if (f <= 0.0 || java.lang.Double.isNaN(f) || java.lang.Double.isNaN(c))
			return 1.0
		
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
		return 1 - chi.cumulativeProbability(DianaUtil.fishers(pvalues) / c)
		//group.rtProb 	= state.rtProb((group.iend + group.istart) / 2)
		//return group
	}
	
	
	
	
	/*
	def getFragmentPvalues():Seq[Double] = {
		var ps = new ArrayBuffer[Double]
		for (di <- 0 until state.numChromatograms)
			ps += 1 / state.targetRatioTable(0)(di).mean
		var sum = ps.sum
		return ps.map(_ / sum)
	}
	
	
	def pvalueForFragments(g:PCGroup, intensities:Seq[Array[Double]]):PCGroup = {
		val p0s	= getFragmentPvalues
		val xs	= intensities.map(_.slice(g.istart, g.iend).sum / 10)
		val n 	= xs.sum
		
		val s = new ArrayBuffer[Double]
		g.fragmentPvalues = p0s.zip(xs).map(t => {
			val p0 = t._1
			val x = t._2
			val my 		= n * p0
			val sigma 	= math.max(math.sqrt(n*p0*(1-p0)), 0.05 * my)
			s += sigma
			math.max(0.001, math.min(0.999, NormalDistribution.cdf((x - my) / sigma)))
		})
		
		g.sigmas = s
		g.xs = xs
		
		val chi = new ChiSquaredDistribution(2*state.numChromatograms)
		g.pEstimates = 1 - chi.cumulativeProbability(DianaUtil.fishers(g.fragmentPvalues))
		
		return g
	}
	*/
	
	
	def corrScore(g:PCGroup, est:GroupEstimation):Double = {
		val xs = est.estimates.map(Filter.savitzkyGolay9 _)//intensities.map(_.slice(g.istart, g.iend))
		val corrs = new ArrayBuffer[Double]
		
		Ratios.iterate(xs.length, (ri, ui, di) => {
			corrs += Stats.pearsonCorrelation(xs(ui), xs(di))
		})
		
		//g.corrs = corrs
		var corrScore = StatUtils.mean(corrs.toArray)
		if (java.lang.Double.isNaN(corrScore)) corrScore = 0.0
		corrScore
	}
}
