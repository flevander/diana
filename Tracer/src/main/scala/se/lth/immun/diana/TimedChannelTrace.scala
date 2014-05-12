package se.lth.immun.diana

class TimedChannelTrace[T](
		val id:T,
		aggr:Int = 1,
		mean:Boolean = false
) {
	val intensity = 
		if (aggr == 1) 	TraceBuilder()
		else			AggrTraceBuilder(aggr, mean)
	val time = 
		if (aggr == 1) 	TraceBuilder()
		else			AggrTraceBuilder(aggr, mean)
		
	override def toString = "TimedChannelTrace["+id+"]"
}