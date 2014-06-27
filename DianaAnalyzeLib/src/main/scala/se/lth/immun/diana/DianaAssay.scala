package se.lth.immun.diana

import se.lth.immun.traml.ghost._

class DianaAssay(
		val pepCompRef:String,
		val pepCompMz:Double,
		val pepCompCharge:Int,
		val expectedRT:Double,
		val transitions:Seq[GhostTransition],
		val targets:Seq[GhostTarget]
) {

}