package se.lth.immun.diana

import akka.actor.Actor
import akka.actor.Props
import akka.actor.PoisonPill
import se.lth.immun.traml.ghost._
import grizzled.slf4j.Logging

class TopActor(params:DianaScorerParams) extends Actor with Logging {
	
	import AnalysisActor._
	import EsvWriteActor._
	
	var timing:Timing = null
    var pcGroups:Option[Seq[Seq[DianaPeptideCandidate]]] = None
    
    val analysisActor = context.actorOf(Props[AnalysisActor], "analysis-actor")
    var analysisTime = 0L
    var traml:GhostTraML = _
    
    def pscore(d:DianaPeptideCandidate):Double = d.fragmentRankPcsRatioProb
    
    def receive = {
		case af:AnalyzeFile =>
			analysisActor ! af
			traml = af.traml
		
		case x:String =>
			print(x)
			
		case AnalysisDone(candidates, t) =>
			println("analysis done received")
			pcGroups = Some(
					candidates.map(_.sortBy(pscore _))
						.sortBy(l => pscore(l.head))
						)
			timing = t
			analysisTime = System.currentTimeMillis - params.startTime - params.chromReadTime - params.tramlReadTime
			if (params.outEsv)
				initiateEsvWrite
			else
				wrapUp
		
		case FileWritten() =>
			wrapUp
		
		case ErrorMsg(msg) =>
			error(msg)
			if (params.verbose)
				println(msg)
		
		case x =>
			println("SOMETHING UNEXPECTED!")
			println(x)
	}
	
	
	
	def initiateEsvWrite = {
		val w = context.actorOf(Props(new EsvWriteActor(params)), "esv-writer")
		w ! Write(traml, pcGroups.get, analysisTime)
	}
	
	
	
	def wrapUp = {
		val after = System.currentTimeMillis
		if (!params.profiling)
			println("  time taken: "+niceTiming(after-params.startTime))
		else {
			println("     total time: "+niceTiming(after-params.startTime))
			println("chrom read time: "+niceTiming(params.chromReadTime))
			println("traml read time: "+niceTiming(params.tramlReadTime))
			println("  analysis time: "+niceTiming(analysisTime))
			println(" out write time: "+niceTiming(params.esvWriteTime))
			val runt = Runtime.getRuntime
			println("      heap size: "+(runt.totalMemory / 1000000) + " / "+(runt.maxMemory / 1000000) + " Mb")
			println()
			println("-- average subtimes per assay --")
			println("           chrom map time: %10.6f ms".format(timing.chromMapTime))
			println("diana ms1 input calc time: %10.6f ms".format(timing.dianaMs1InputTime))
			println("diana ms2 input calc time: %10.6f ms".format(timing.dianaMs2InputTime))
			println("      diana analysis time: %10.6f ms".format(timing.dianaAnalysisTime))
		}
		
		context.system.shutdown
	}
	
	
	def niceTiming(t:Long):String = {
		
		val ms = t % 1000
		var x = t / 1000
		val s = x % 60
		x = x / 60
		val m = x & 60
		x = x / 60
		val h = x % 24
		val d = x / 24
		val str = "%d days %02d:%02d:%02d.%03d".format(d, h, m, s, ms)
		val init = str.takeWhile(c => !c.isDigit || c == '0')
		init.replace('0', '_') + str.drop(init.length)
	}
	
}