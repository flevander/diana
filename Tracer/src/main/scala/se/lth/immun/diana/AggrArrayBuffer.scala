package se.lth.immun.diana

import collection.mutable.ArrayBuffer

class AggrArrayBuffer[T](
		val aggr: Int, 
		val mean:Boolean = false
) extends ArrayBuffer[Double] {
	
	private var buffSum = 0.0
	private var buffCount = 0
	
	override def +=(elem: Double): this.type = {
		buffSum += elem
		buffCount += 1
		if (buffCount == aggr) {
			 super.+=(
					 if (mean) buffSum / aggr
					 else buffSum
					 )
			buffCount = 0
			buffSum = 0
		}
		this
	}
}