package se.lth.immun.diana

import scala.collection.immutable.TreeMap
import scala.collection.mutable.ArrayBuffer

object ChannelLookup {
	case class Range(start:Double, end:Double) {
		override def toString = start +"-"+end
	}
}

class ChannelLookup[C](
		channels:Seq[C],
		channelPos:C => Double
) {
	import ChannelLookup._
	
	val channelMap = {
		var map = TreeMap[Double, ArrayBuffer[C]]() ++ channels.map(channelPos).toSet.map((d:Double) => (d, new ArrayBuffer[C]))
		for (c <- channels) 
			map(channelPos(c)) += c 
		map
	}
	
	
	def fromRanges(ranges:Seq[Range]) = 
		ranges.map(r => channelMap.range(r.start, r.end)).fold(Nil.toMap)(_ ++ _).values.flatten
}