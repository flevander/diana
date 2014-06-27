package se.lth.immun.diana

import se.lth.immun.mzml.ghost._
import se.lth.immun.math.Ratios

import scala.collection.mutable.ArrayBuffer

object DianaInput {
	
	case class Channel(mz:Double, id:String)
	case class Ratio(id1:Int, id2:Int, ratio:Double)
	
	def inverseRatio(r:Ratio) = Ratio(r.id2, r.id1, 1 / r.ratio)
	val IDENTITY_RATIO = Ratio(0, 0, 1)
	
	/**
	 * This function assumes matching chroms and trans sequences, and that no chromatograms are missing.
	 */
	def fromTransitions(chroms:Seq[XChromatogram], assay:DianaAssay, times:Seq[Double] = Nil):DianaInput = {
		var refRatios = new ArrayBuffer[Ratio]
	
	    Ratios.iterate(assay.transitions.length, 
	    			(ri, i, j) => {
	    				refRatios += Ratio(i, j, assay.transitions(i).intensity / assay.transitions(j).intensity) 
	    			}
	    		)
	    
	    val cg = new XChromatogramGroup(assay.transitions.head.q1)
		cg.chromatograms ++= chroms
	    		
	    new DianaInput(
	    		assay.pepCompRef,
	    		assay.pepCompMz,
	    		assay.pepCompCharge,
	    		cg.resample(times),
	    		assay.transitions.map(t => Channel(t.q1, t.id)),
	    		refRatios)
	}
	
	
	
	
	/**
	 * This function assumes matching chroms and trans sequences, and that no chromatograms are missing.
	 */
	def fromTargets(chroms:Seq[XChromatogram], assay:DianaAssay, times:Seq[Double] = Nil):DianaInput = {
		var refRatios = new ArrayBuffer[Ratio]
	
	    Ratios.iterate(assay.targets.length, 
	    			(ri, i, j) => {
	    				refRatios += Ratio(i, j, assay.targets(i).intensity / assay.targets(j).intensity) 
	    			}
	    		)
	    
	    val cg = new XChromatogramGroup(assay.transitions.head.q1)
		cg.chromatograms ++= chroms
		
	    new DianaInput(
	    		assay.pepCompRef,
	    		assay.pepCompMz,
	    		assay.pepCompCharge,
	    		cg.resample(times),
	    		assay.targets.map(t => Channel(t.q1, t.id)),
	    		refRatios)
	}
}

class DianaInput(
		val ref:String,
		val mz:Double,
		val charge:Int,
		val cg:XChromatogramGroup,
		val channels:Seq[DianaInput.Channel],
		val refRatios:Seq[DianaInput.Ratio]
		
) {

}