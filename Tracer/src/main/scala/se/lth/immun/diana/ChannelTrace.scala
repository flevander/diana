package se.lth.immun.diana

import se.lth.immun.collection.numpress.NumSlofArray

class ChannelTrace[T](
		val id:T,
		aggr:Int = 1,
		mean:Boolean = false
) {
	val intensity = 
		new AggrBuilder(new NumSlofArray(22, 20), aggr, mean)
}