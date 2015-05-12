package se.lth.immun.diana

import collection.mutable.HashMap
import collection.mutable.HashSet
import collection.mutable.Queue
import collection.mutable.ArrayBuffer
import util.Random

import akka.actor.Actor
import akka.actor.ActorRef
import akka.actor.Props

import java.io.File
import se.lth.immun.mzml.ghost._
import se.lth.immun.traml.ghost._

import se.jt.CLIBar


object AnalysisActor {
	
	case class AnalyzeFile(
			traml:GhostTraML,
			params:DianaScorerParams
			)
			
	case class AnalysisDone(
			candidates:Seq[Seq[DianaPeptideCandidate]],
			t:Timing
			)
		
	case class Timing(
			chromMapTime:Double, 
			dianaMs1InputTime:Double, 
			dianaMs2InputTime:Double, 
			dianaAnalysisTime:Double)
}

class AnalysisActor extends Actor {

	import AnalysisActor._
	import AssayAnalyzer._
	
	var swathGmzML:GhostMzML	= _
	var cliBar		= new CLIBar(30, true)
	
	var t1 	= 0L
	var t2 	= 0L
	var t3 	= 0L
	var t4 	= 0L
	
	val assayQueue = new Queue[DianaAssay]
	var runningAssays = new HashSet[DianaAssay]()
	val results = new ArrayBuffer[Seq[DianaPeptideCandidate]]()
	var reportAssays:Set[DianaAssay] = _
	
	var customer:ActorRef = _
	var nAssays = 0
	var aCount = 0
	var params:DianaScorerParams = _
	
	def receive = {
		case AnalyzeFile(traml, p) =>
			customer = sender
			params = p
			analyzeFile(p.swathFile, traml)
			
		case AssayAnalyzed(da, res, timing) =>
			results += res
			
			t1 += timing.t1
			t2 += timing.t2
			t3 += timing.t3
			t4 += timing.t4
						
			aCount += 1
			if (!params.verbose)
				customer ! cliBar.update(aCount, nAssays)
			/*else
				customer ! "%4d ms\t+%d\t%s\n".format(timing.t4, 
								da.pepCompCharge, da.pepCompRef)
			*/
			runningAssays -= da
			if (assayQueue.nonEmpty) 
				startAssayAnalysis(assayQueue.dequeue)
			
			if (runningAssays.isEmpty) 
				wrapUpAnalysis
			
			context.stop(sender)
		}
	
	def analyzeFile(
			mzmlFile:File,
			traml:GhostTraML
	):Unit = {
		val (reader, binFile, binFileChannel) = GhostMzML.getReaders(mzmlFile)
		sender ! " reading chromatogram mzML file...\n"
		val before = System.currentTimeMillis
		swathGmzML = GhostMzML.fromFile(reader, true, binFileChannel)
		params.chromReadTime = System.currentTimeMillis - before
		
		val assays = getAssays(traml)
		if (params.reportSeed >= 0)
		Random.setSeed(params.reportSeed)
		reportAssays = Random.shuffle(assays).take(params.nReport).toSet
		
		if (params.nReport > 0) {
			val qcDir = new File("qc")
			if (!qcDir.exists)
				qcDir.mkdir
		}
		
		aCount 	= 0
		nAssays = assays.length
		
		cliBar.reset
		if (!params.verbose) {
			sender ! cliBar.reference + "\n"
			sender ! cliBar.update(0)
		}
		
		assayQueue.enqueue(assays.drop(params.concurrency):_*)
		for (a <- assays.take(params.concurrency)) 
			startAssayAnalysis(a)
	}
	
	
	
	def startAssayAnalysis(a:DianaAssay) = {
		import AssayAnalyzer._
		
		runningAssays += a
		val actor = context.actorOf(Props(new AssayAnalyzer(params)))
		
		//val transChroms = a.transitions.map(t => swathGmzML.chromatograms.find(_.id == t.id))
		//val targChroms = a.targets.map(t => swathGmzML.chromatograms.find(_.id == t.id))
		actor ! AnalyzeAssay(a, swathGmzML, reportAssays.contains(a))//transChroms, targChroms)
		if (params.verbose)
			customer ! "+%d %s\n".format(a.pepCompCharge, a.pepCompRef)
		
	}
	
	
	
	def wrapUpAnalysis = {
		val n = nAssays.toDouble
		customer ! AnalysisDone(results, 
				AnalysisActor.Timing(t1 / n, t2 / n, t3 / n, t4 / n))
	}
	
	
	
	def getAssays(traml:GhostTraML):Seq[DianaAssay] = {
		
		val pcBuff = new HashMap[(String, Int), (Double, ArrayBuffer[GhostTransition], ArrayBuffer[GhostTarget])]
		for (t <- traml.transitions) {
			val z = 
				if (t.q1z != 0) t.q1z
				else if (traml.peptides.contains(t.peptideRef)) {
					import PeptideParser._
					val seq = traml.peptides(t.peptideRef).sequence
					PeptideParser.parseSequence(seq) match {
						case UniModPeptide(p) => math.round(p.monoisotopicMass() / t.q1).toInt
						case Unparsable(seq) => 
							val errMsg = "WARN: could parse peptide '%s' for mass calculation".format(seq)
							customer ! ErrorMsg(errMsg)
							0
					}
				} else {
					val errMsg = "WARN: no charge given for compound '%s' with mz %.3f".format(t.compoundRef, t.q1)
					customer ! ErrorMsg(errMsg)
					0
				}
			val key = (if (t.peptideRef != null) t.peptideRef else t.compoundRef, z)
			if (!pcBuff.contains(key))
				pcBuff += key -> (t.q1, new ArrayBuffer, new ArrayBuffer)
			pcBuff(key)._2 += t
		}
		for (t <- traml.includes) {
			val z = 
				if (t.q1z != 0) t.q1z
				else if (traml.peptides.contains(t.peptideRef)) {
					import PeptideParser._
					val seq = traml.peptides(t.peptideRef).sequence
					PeptideParser.parseSequence(seq) match {
						case UniModPeptide(p) => math.round(p.monoisotopicMass() / t.q1).toInt
						case Unparsable(seq) => 
							val errMsg = "WARN: could parse peptide '%s' for mass calculation".format(seq)
							sender ! ErrorMsg(errMsg)
							0
					}
				} else {
					val errMsg = "WARN: no charge given for compound '%s' with mz %.3f".format(t.compoundRef, t.q1)
					sender ! ErrorMsg(errMsg)
					0
				}
			val key = (if (t.peptideRef != null) t.peptideRef else t.compoundRef, z)
			if (!pcBuff.contains(key))
				pcBuff += key -> (t.q1, new ArrayBuffer, new ArrayBuffer)
			pcBuff(key)._3 += t
		}
		
		pcBuff.toSeq.map(kv => {
			val ((ref, z), (tmz, trans, targs)) = kv
			
			
			val expectedRT = trans.map(t => (t.rtEnd + t.rtStart)/2).sum / trans.length
			new DianaAssay(ref, tmz, z, expectedRT, trans, targs)
		})
	}
}