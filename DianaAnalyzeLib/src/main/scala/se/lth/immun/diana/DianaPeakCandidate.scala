package se.lth.immun.diana

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.distribution.ChiSquaredDistribution

import se.lth.immun.math.Matrix
import se.lth.immun.markov.MarkovChainDistribution
import DianaInput._

object DianaPeakCandidate {
	
	
	class PC(
			val istart:Int,
			val iapex:Int,
			val iend:Int
	) {
		var icurve = -1
	}
	
	
	
	
	class Validation(
		var icurve:Int,
		width:Int
	) {
		var ok 		= new Array[Boolean](width)
		var flag 	= false
		for (i <- 0 until width) ok(i) = false
	}
	
	
	
	
	
	object GroupValidation {
		
		def apply(nValids:Int, iwidth:Int) = {
			var groupValids = new GroupValidation(
				new Array[Validation](nValids),
				Matrix.get2d[Validation](nValids).map(_.map(v => new Validation(0, iwidth)))
			)
			for (i <- 0 until nValids)
				groupValids.valids(i) = new Validation(i, iwidth)
			groupValids
		}
	}
	
	class GroupValidation(
		var valids:Array[Validation],
		var rValids:Array[Array[Validation]]
	) {}
	
	
	
	
	
	object GroupEstimation {
		
		def apply(nChrom:Int, niwidth:Int) = {
			val est = new GroupEstimation
			
			est.estimateApex = 0.0
			est.estimates 	= new Array[Array[Double]](nChrom)
	        for (e <- 0 until nChrom) 
	        	est.estimates(e) = new Array[Double](niwidth)
	        
	        est.areas = new Array(nChrom)
	        est
		}
	}
	
	class GroupEstimation {
		var estimates:Array[Array[Double]] = null
		var estimateApex:Double = -1.0
		var iEstimateApex:Int = -1
		var areas:Array[Double] = null
			
		def area = if (areas == null || areas.isEmpty) 0.0 else areas.sum
	}
	
	
	
	
	
	class PCGroup {
		var pcs 	= new ArrayBuffer[PC]
		//var rtProb 		= 1.0
		//var ratioProb 	= 1.0
		//var pEstimates	= 1.0
		//var fragmentPvalues:Seq[Double] = null
		//var sigmas:Seq[Double] = null
		//var xs:Seq[Double] = null
		//var corrScore:Double = 0.0
		//var corrs:Seq[Double] = null
		//var pvalue 	= -1.0
		//var qvalue 	= -1.0
		var pvalueTable:Array[Array[Double]] 	= null
		//var validations:Array[Validation]		= null
		//var ratioValidations:Array[Array[Validation]] = null
		
		override def toString = istart + "-"+ iend
	
		def iend 	= StatUtils.percentile(pcs.map(_.iend.toDouble).toArray, 50).toInt //pcs.map(_.iend).max
		def istart 	= StatUtils.percentile(pcs.map(_.istart.toDouble).toArray, 50).toInt //pcs.map(_.istart).min
		def iwidth 	= iend - istart
		
		var nistart = 0
		var niend = 0
		def niwidth = niend - nistart
		
		def getQuota(ui:Int, di:Int, gv:GroupValidation) = {
			val ok = gv.rValids(ui)(di).ok
			val count = 
				math.max(
						0, 
						ok.slice(nistart, niend).count(b => b)-1
					)
			count.toDouble / niwidth
		}
		
		
		
		def estimateAndIntegrate(
				state:DianaChromatogramState, 
				y:Array[Array[Double]],
				gv:GroupValidation
		):GroupEstimation = {
	        def stable(r:Ratio):Boolean = true//0.5 * r.mean > r.stdDev
	
	        val L		= state.numChromatograms
	        val t0		= istart
	        val est 	= GroupEstimation(L, niwidth)
	        
	
	        for (t <- nistart until niend) {
	            for (b <- 0 until L) {
	                try {
	                	if (gv.valids(b).ok(t))
		                    est.estimates(b)(t - nistart) = y(b)(t0 + t)
		                else {
		                    var localEstimates = List[Double]()
		                    for (e <- 0 until L) {
		                        if (
		                        			e != b 
		                        		&& 	gv.valids(e).ok(t) 
		                        		&& 	stable(state.targetRatioTable(b)(e))
		                        ) {
		                        	var multiplier = state.targetRatioTable(b)(e).ratio
		                        	localEstimates = localEstimates :+ (y(e)(t0 + t) * multiplier)
		                        }
		                    }
		                    est.estimates(b)(t - nistart) = 
		                    	if (localEstimates.isEmpty)	 y(b)(t0 + t) 
		                    	else 
		                    		math.min(StatUtils.mean(localEstimates.toArray), y(b)(t0 + t))
		                }
	                } catch {
	                	case e:Exception => e.printStackTrace()
	                }
	            }
	        }
	        
	        
	        for (ifrag <- 0 until L) {
	            est.areas(ifrag) = 0
	            var vals = y(ifrag)
	            for (et <- 0 until niwidth) {
	            	val x = 
		                if (vals(t0 + nistart + et) <= est.estimates(ifrag)(et))
		                	vals(t0 + nistart + et)
						else 
							est.estimates(ifrag)(et)
					est.areas(ifrag) += x
					if (x > est.estimateApex) {
						est.estimateApex = x
						est.iEstimateApex = t0 + nistart + et
					}
	            }
	        }
	        
	        est
		}
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
		val bs		= if (_binSize > 0) _binSize else DianaUtil.DEFAULT_BIN_SIZE
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
			pcs:Seq[Seq[PC]],
			findEdges:Seq[PC] => Array[List[Int]]
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
}
