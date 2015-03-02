package se.lth.immun.diana

import se.lth.immun.collection.numpress.NumLinArray
import se.lth.immun.collection.numpress.NumSlofArray

class ChannelGroupTrace[G,SG](
		val id:G,
		val subChannelIds:Seq[SG]
) {
	
	val time = 
		new NumLinArray(10000.0, 20)

	val subChannels =
		subChannelIds.map(id =>
			id -> new NumSlofArray(22.0, 20)
		).toMap
}