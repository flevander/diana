package se.lth.immun.diana

import scala.collection.mutable.ArrayBuilder
import scala.collection.mutable.HashMap
import se.lth.immun.traml.ghost.GhostTransition

object TransitionPartitioning {
	case class Q1Window(min:Double, max:Double) {
		def has(q1:Double) = q1 >= min && q1 < max
	}
	
	class PartitionExtract(val chroms:Array[ChromExtract[GhostTransition]]) {
		val times = ChromBuilder()
	}
}
class TransitionPartitioning(
		val transitions:Seq[ChromExtract[GhostTransition]]
) {
	import TransitionPartitioning._
	
	val parts = HashMap[Q1Window, PartitionExtract]()
	
	def getTransitions(p:Q1Window) = {
		if (!(parts contains p))
			parts += p -> new PartitionExtract(transitions.filter(ce => p.has(ce.id.q1)).toArray)
		parts(p)
	}
}