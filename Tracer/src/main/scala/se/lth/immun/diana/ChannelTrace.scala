package se.lth.immun.diana

import se.lth.immun.collection.numpress.NumSlofArray
import collection.mutable.ArrayBuffer

class ChannelTrace[T](
		val id:T,
		aggr:Int = 1,
		mean:Boolean = false,
		intFP:Double = 22.0
) {
	val intensity = new AggrBuilder(new NumSlofArray(intFP, 20), aggr, mean)
}