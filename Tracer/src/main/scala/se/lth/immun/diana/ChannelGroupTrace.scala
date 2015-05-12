package se.lth.immun.diana

import se.lth.immun.collection.numpress.NumLinArray
import se.lth.immun.collection.numpress.NumSlofArray

import collection.mutable.ArrayBuffer

class ChannelGroupTrace[G,SG](
		val id:G,
		val subChannelIds:Seq[SG],
		timeFP:Double = 100000.0,
		intFP:Double = 22.0
) {
	
	val time = new NumLinArray(timeFP, 20)

	val subChannels =
		subChannelIds.map(id =>
			id -> new NumSlofArray(intFP, 20)
		).toMap
}