package se.lth.immun.diana

import scala.collection.mutable.ArrayBuilder
import scala.collection.mutable.WrappedArray

object AggrTraceBuilder {
	def apply(aggr:Int, mean:Boolean) = new AggrTraceBuilder(aggr, mean)
}

class AggrTraceBuilder(val aggr: Int, val mean:Boolean = false, sizeIncrement: Double = 2.0) extends TraceBuilder(sizeIncrement) {

	private var buffSum = 0.0
	private var buffCount = 0

	override def +=(elem: Double): this.type = {
		buffSum += elem
		buffCount += 1
		if (buffCount == aggr) {
			ensureSize(size + 1)
			elems(size) = 
				if (mean) buffSum / aggr
				else buffSum
			size += 1
			buffCount = 0
			buffSum = 0
		}
		this
	}

	override def ++=(xs: TraversableOnce[Double]): this.type = xs match {
		case xs: WrappedArray.ofDouble =>
			ensureSize(this.size + xs.length)
			for (x <- xs) {
				buffSum += x
				buffCount += 1
				if (buffCount == aggr) {
					elems(size) = 
						if (mean) buffSum / aggr
						else buffSum
					size += 1
					buffCount = 0
					buffSum = 0
				}
			}
			this
		case _ =>
			super.++=(xs)
	}
	
	override def toString = "AggrTraceBuilder"
}