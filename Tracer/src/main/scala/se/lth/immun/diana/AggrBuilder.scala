package se.lth.immun.diana

import scala.collection.mutable.Builder
import se.lth.immun.collection.numpress.NumpressArray

class AggrBuilder(
		val a:NumpressArray, 
		val aggr: Int, 
		val mean:Boolean = false
) {

	private var buffSum = 0.0
	private var buffCount = 0
	
	def +=(elem: Double): this.type = {
		buffSum += elem
		buffCount += 1
		if (buffCount == aggr) {
			 a += 	(
					 if (mean) buffSum / aggr
					 else buffSum
					 )
			buffCount = 0
			buffSum = 0
		}
		this
	}
}