package se.lth.immun.diana

import akka.actor.Actor
import akka.actor.Props
import se.lth.immun.mzml.ghost._

import se.lth.immun.math.Ratios
import se.lth.immun.signal.Filter

import se.lth.immun.mzml.ghost.XChromatogram

object AssayAnalyzer {
	
	case class Timing(t1:Long, t2:Long, t3:Long, t4:Long)
	case class AnalyzeAssay(
			da:DianaAssay, gmzML:GhostMzML, writeReport:Boolean)
	//		transChroms:Seq[Option[GhostChromatogram]], 
	//		targChroms:Seq[Option[GhostChromatogram]]
	//		)
	case class AssayAnalyzed(da:DianaAssay, res:Seq[DianaPeptideCandidate], t:Timing)
	
	val signalP = DianaSignalProcessor.getDefault
	val carrierAreaCutoff = 25.0
}




class AssayAnalyzer(params:DianaScorerParams) extends Actor {

	import AssayAnalyzer._
	
	var lastTime = 0L
	def t0 = lastTime = System.currentTimeMillis
	def dt = {
		val t = System.currentTimeMillis
		val dt = t - lastTime
		lastTime = t
		dt
	}
	
	
	def receive = {
		case AnalyzeAssay(da, gmzML, writeReport) =>
			
			t0
			val transChroms = da.transitions.map(t => gmzML.chromatograms.find(c => c.id == t.id && c.intensities.length > 0))
			val targChroms = da.targets.map(t => gmzML.chromatograms.find(c => c.id == t.id && c.intensities.length > 0))
			
			val t1 = dt
			
			if (transChroms.count(_.isDefined) < 2) {
				val errMsg = "ERROR during analysis of '%s %.3f'".format(da.pepCompRef, da.pepCompMz)
        		sender ! ErrorMsg(errMsg) 
            	sender ! AssayAnalyzed(da, List(new DianaPeptideCandidate(da)), Timing(t1, 0, 0, 0))
			} else {
				
				val dIn = DianaInput.fromTransitions(transChroms.map(_.get.toXChromatogram), da)
				
				val t2 = dt
				
				val iDIn =
					if (targChroms.nonEmpty)
						Some(DianaInput.fromTargets(
								targChroms.map(_.get.toXChromatogram), 
								da, dIn.cg.chromatograms.head.times
							))
					else None
					
				val t3 = dt
				
				val carriers = analyzeCG(dIn, iDIn, da)
				val t4 = dt
				
				val res =
					if (carriers.isEmpty) List(new DianaPeptideCandidate(da))
			        else carriers.map(c => 
			        	DianaPeptideCandidate(da, dIn.cg, iDIn.map(_.cg), 
			        			c, carriers.length))
			    sender ! AssayAnalyzed(da, res, Timing(t1, t2, t3, t4))
			    
			    if (writeReport) {
			    	import AssayReportActor._
		
					val reporter = context.actorOf(Props(new AssayReportActor(params)))
					reporter ! WriteReport(da, carriers, dIn, iDIn)
			    }
			}
		}
	
	
	
	
	
	val fragmentState 		= new DianaChromatogramState(params.relRatioUpperBound, 
													Double.NaN, Double.NaN, Double.NaN, 
													params.signalBinSize)
	val fragmentEvaluator 	= new DianaPCEvaluator(fragmentState)
	var findEdges			= fragmentEvaluator.findEdges _
	
	val isotopeState 		= new DianaChromatogramState(1.5)
	val isotopeEvaluator 	= new DianaPCEvaluator(isotopeState)
	
	def analyzeCG(dIn:DianaInput, iDIn:Option[DianaInput], assay:DianaAssay):Seq[Carrier] = {
		
		val fChroms = dIn.cg.chromatograms
		val iChroms:Seq[XChromatogram] = iDIn match {
			case Some(x) => x.cg.chromatograms 
			case None => Nil
		}
		
		val allChroms 	= fChroms ++ iChroms
		if (allChroms.isEmpty) return Nil 
			
		def getPCs(y:DianaSignalProcessor.SmoothAndBase) = {	
			var dy 		= signalP.getDerivate(y.smooth)
			var ddy		= signalP.getDerivate(dy)
			
			DianaPeakCandidate.findPCs(y.smooth, dy, ddy, y.base)
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
		val carriers	= DianaPeakCandidate.groupPCs(pcs, findEdges)
							.map(g => new Carrier(assay, g, fChroms.length, iChroms.length))
							.filter(nonRetarded)
		
		
		val filteredCarriers = calculateFragmentScores(dIn, smooths, assay, carriers)
		if (iDIn.nonEmpty)
			calculateIsotopeScores(iDIn.get, assay, filteredCarriers)
		
		for (c <- filteredCarriers)
			c.rtProb = fragmentState.rtProb(c.fragmentEstimation.iEstimateApex)
		
		return filteredCarriers
	}
	
	
	
	
	def calculateFragmentScores(
			dIn:DianaInput, 
			smooths:Seq[DianaSignalProcessor.SmoothAndBase],
			assay:DianaAssay, 
			carriers:Seq[Carrier]
	) = {
		val chroms 	= dIn.cg.chromatograms
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
				ratios.getTargetRatioTable(DianaInput.IDENTITY_RATIO)(dIn.refRatios.toArray, DianaInput.inverseRatio),
				assay.expectedRT,
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
		
		val ret = carriers.filter(_.isRemotelyUnlikely(params.pCutoff))
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
			dIn:DianaInput, 
			assay:DianaAssay, 
			carriers:Seq[Carrier]
	) = {
		val chroms 	= dIn.cg.chromatograms
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
				ratios.getTargetRatioTable(DianaInput.IDENTITY_RATIO)(dIn.refRatios.toArray, DianaInput.inverseRatio),
				assay.expectedRT,
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