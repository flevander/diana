package se.lth.immun.diana

import se.lth.immun.anubis.Ratio
import se.lth.immun.anubis.AnubisInput
import se.lth.immun.anubis.ReferencePrecursor
import se.lth.immun.math.Ratios
import se.lth.immun.signal.Filter

import se.lth.immun.mzml.ghost.XChromatogram

object AnalyzeCG {
	
	import SwathPeakCandidate._
	import SwathSignalProcessor._


	var pValueCutoff 	= 0.99
	val signalP 		= SwathSignalProcessor.getDefault
	
	val fragmentState 		= new SwathChromatogramState(1.5)
	val fragmentEvaluator 	= new SwathPCEvaluator(fragmentState)
	var findEdges			= fragmentEvaluator.findEdges _
	
	val isotopeState 		= new SwathChromatogramState(1.5)
	val isotopeEvaluator 	= new SwathPCEvaluator(isotopeState)
	
	val carrierAreaCutoff 	= 25.0

	
	
	def apply(aIn:AnubisInput, iAIn:Option[AnubisInput], pc:ReferencePrecursor):Seq[Carrier] = {
		
		val before 		= System.currentTimeMillis
		
		val fChroms = aIn.cg.chromatograms
		val iChroms:Seq[XChromatogram] = iAIn match {
			case Some(x) => x.cg.chromatograms 
			case None => Nil
		}
		
		val allChroms 	= fChroms ++ iChroms
		if (allChroms.isEmpty) return Nil 
			
		def getPCs(y:SwathSignalProcessor.SmoothAndBase) = {	
			var dy 		= signalP.getDerivate(y.smooth)
			var ddy		= signalP.getDerivate(dy)
			
			SwathPeakCandidate.findPCs(y.smooth, dy, ddy, y.base)
		}
		
		def nonRetarded(c:Carrier) = {
			val nFrags = c.fragmentPcs.length 
			nFrags > 1 || 
			nFrags > 0 && {
				val pc = c.fragmentPcs.head
				fChroms(pc.icurve).intensities.slice(pc.istart, pc.iend).sum > carrierAreaCutoff
			}
		}
		val smooths		= allChroms.map(xc => signalP.getSmoothAndBase(xc.intensities.toArray))
		val pcs 		= smooths.map(xc => getPCs(xc))
		val carriers	= SwathPeakCandidate.groupPCs(pcs, findEdges)
							.map(g => new Carrier(pc, g, fChroms.length, iChroms.length))
							.filter(nonRetarded)
		
		
		val filteredCarriers = calculateFragmentScores(aIn, smooths, pc, carriers)
		if (iAIn.nonEmpty)
			calculateIsotopeScores(iAIn.get, pc, filteredCarriers)
		
		for (c <- filteredCarriers)
			c.rtProb = fragmentState.rtProb(c.fragmentEstimation.iEstimateApex)
		
		if (!DianaScorer.quiet)
			println("%7.2f\t%4d ms\t%s".format(pc.mz, 
				System.currentTimeMillis-before, pc.peptideSequence))
		
		return filteredCarriers
	}
	
	
	
	def calculateFragmentScores(
			aIn:AnubisInput, 
			smooths:Seq[SmoothAndBase],
			pc:ReferencePrecursor, 
			carriers:Seq[Carrier]
	) = {
		val chroms 	= aIn.cg.chromatograms
		val l 		= chroms.length
		val t 		= chroms(0).times.toArray
		
		var reduced = Filter.baseLineReduce(
        					chroms.toArray.map(tr => tr.intensities.toArray)
        				)
        var filtered 	= reduced.map(ds => Filter.savitzkyGolay9(ds))
        var ratios 		= new Ratios(filtered)
		fragmentState.setChromatogram(
				l, 
				ratios, 
				ratios.getTargetRatioTable(Ratio.ONE)(aIn.refRatios, Ratio.INVERSE),
				pc.retentionTime.peak,
				t
			)
		
		for (c <- carriers)
			c.fragmentValids = fragmentEvaluator.validateGroup(c.g, c.fragmentPcs, smooths)
		
		if (carriers.nonEmpty)
			fragmentState.calculateChromatogramStatistics(carriers.map(c => (c.g, c.fragmentValids)))
		
		for (c <- carriers) {
			c.fragmentRankAllRatioProb = fragmentEvaluator.nullRatioProb(
									c.g, 
									c.fragmentValids,
									fragmentState.statsAll,
									"rank")
			c.fragmentRankPcsRatioProb = fragmentEvaluator.nullRatioProb(
									c.g, 
									c.fragmentValids,
									fragmentState.statsPcs,
									"rank")
			c.fragmentMarkovAllRatioProb = fragmentEvaluator.nullRatioProb(
									c.g, 
									c.fragmentValids,
									fragmentState.statsAll,
									"markov")
			c.fragmentMarkovPcsRatioProb = fragmentEvaluator.nullRatioProb(
									c.g, 
									c.fragmentValids,
									fragmentState.statsPcs,
									"markov")
		}
		
		val ret = carriers.filter(_.isRemotelyUnlikely(pValueCutoff))
		for (c <- ret) {
			c.fragmentEstimation = c.g.estimateAndIntegrate(
									fragmentState, 
									chroms.map(_.intensities.toArray).toArray, 
									c.fragmentValids)
			c.fragmentCorrScore = fragmentEvaluator.corrScore(c.g, c.fragmentEstimation)
		}
		ret
	}
	
	
	
	def calculateIsotopeScores(
			aIn:AnubisInput, 
			pc:ReferencePrecursor, 
			carriers:Seq[Carrier]
	) = {
		val chroms 	= aIn.cg.chromatograms
		val l 		= chroms.length
		val t 		= chroms(0).times.toArray
		val smooths	= chroms.map(xc => signalP.getSmoothAndBase(xc.intensities.toArray))
		
		var reduced = Filter.baseLineReduce(
        					chroms.toArray.map(tr => tr.intensities.toArray)
        				)
        var filtered 	= reduced.map(ds => Filter.savitzkyGolay9(ds))
        var ratios 		= new Ratios(filtered)
		isotopeState.setChromatogram(
				l, 
				ratios, 
				ratios.getTargetRatioTable(Ratio.ONE)(aIn.refRatios, Ratio.INVERSE),
				pc.retentionTime.peak,
				t
			)
		
		for (c <- carriers)
			c.isotopeValids = isotopeEvaluator.validateGroup(c.g, c.isotopePcs, smooths)
		
		if (carriers.nonEmpty)
			isotopeState.calculateChromatogramStatistics(carriers.map(c => (c.g, c.isotopeValids)))
		
		for (c <- carriers) {
			c.isotopeRankAllRatioProb = isotopeEvaluator.nullRatioProb(
									c.g, 
									c.isotopeValids,
									isotopeState.statsAll,
									"rank")
			c.isotopeRankPcsRatioProb = isotopeEvaluator.nullRatioProb(
									c.g, 
									c.isotopeValids,
									isotopeState.statsPcs,
									"rank")
			c.isotopeMarkovAllRatioProb = isotopeEvaluator.nullRatioProb(
									c.g, 
									c.isotopeValids,
									isotopeState.statsAll,
									"markov")
			c.isotopeMarkovPcsRatioProb = isotopeEvaluator.nullRatioProb(
									c.g, 
									c.isotopeValids,
									isotopeState.statsPcs,
									"markov")
			c.isotopeEstimation = c.g.estimateAndIntegrate(
									isotopeState, 
									chroms.map(_.intensities.toArray).toArray, 
									c.isotopeValids)
			c.isotopeCorrScore = isotopeEvaluator.corrScore(c.g, c.isotopeEstimation)
		}
	}
	
	
	
	/*
	def looseFindEdges(
			matchReq:Int,
			state:SwathChromatogramState
		)(
			pcs:Seq[SwathPeakCandidate.PC]
	):Array[List[Int]] = {
		var L 		= pcs.length
		var edges 	= new Array[List[Int]](L)
		val ub 		= state.upperBound
		val lb 		= state.lowerBound
		
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
            	var r 		= state.ratios.getRatio(ui, di)
            	var target 	= state.targetRatioTable(ui)(di).mean
            	var count 	= 0
            	for (ir <- ki.istart until li.iend+1) {
            		if (		r(ir) < target * ub
            				&&	r(ir) > target * lb )
            			count += 1
            	}
                if (count > matchReq)
                	edges(i) = edges(i) :+ k
                
                k += 1
            }
        }
		return edges
	}
	*/
}
