package se.lth.immun.diana

import scala.collection.mutable.ArrayBuilder

class ChromExtract[T](
		val id:T,
		aggr:Int = 1,
		mean:Boolean = false
) {
	val chrom = 
		if (aggr == 1) 	ChromBuilder()
		else			AggrChromBuilder(aggr, mean)
}