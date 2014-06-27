package se.lth.immun.diana

import se.lth.immun.signal.Wavelet
import se.lth.immun.signal.MODWT
import se.lth.immun.signal.WaveletFilter
import se.lth.immun.signal.ITransform
import se.lth.immun.signal.DWaveletLevel

import org.apache.commons.math3.stat.StatUtils

object DianaSignalProcessor {
	
	val MIN_BASELINE = 1.0
	
	def getDefault = 
		new DianaSignalProcessor(
				new WaveletFilter(WaveletFilter.Type.LA8),
				new MODWT,
				Wavelet.Boundary.PERIODIC,
				2,
				20
		)
	
	class SmoothAndBase(
			val smooth:Array[Double],
			val base:Seq[Double]
	) {}
}

class DianaSignalProcessor(
	val filter:WaveletFilter,
	val modwt:ITransform,
	val boundary:Wavelet.Boundary,
	val wlevels:Int,
	val noBaseLine:Int
) {
	import DianaSignalProcessor._
	
	def getSmoothAndBase(x:Array[Double]):SmoothAndBase = {
		var wavelets:Array[DWaveletLevel] = Wavelet.decompose(x, x.length, wlevels, 
											filter, modwt, boundary, null)
		val binSize 		= x.length / noBaseLine
		var smoothLevels 	= Wavelet.smooths(wavelets, wlevels, filter, modwt, boundary, null)
		var y				= smoothLevels(wlevels-1)
		var baseline 		= y.sliding(binSize).map(w => 
									math.max(MIN_BASELINE, StatUtils.percentile(w, 50))
								).toSeq
		new SmoothAndBase(y, baseline)
	}
	
	/**
	 * Computes the derivate of a uniformly sampled vector x of size n
	 * 
	 * @return 	Array of size n-1 with derivates for the midpoints for 
	 * 			the point halfway between the original samplings 
	 */
	def getDerivate(x:Array[Double]):Array[Double] = 
		return x.zip(x.tail).map(tu => tu._2 - tu._1)
}
