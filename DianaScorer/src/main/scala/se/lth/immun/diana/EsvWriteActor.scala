package se.lth.immun.diana

import akka.actor.Actor

import se.lth.immun.esv._
import se.lth.immun.traml.ghost._

import java.io.File
import java.io.FileWriter
import java.io.BufferedWriter
import java.util.Locale

object EsvWriteActor {
	
	case class Write(
			traml:GhostTraML,
			pcGroups:Seq[Seq[DianaPeptideCandidate]], 
			analysisTime:Long)
	case class FileWritten()
}


class EsvWriteActor(params:DianaScorerParams) extends Actor {

	import EsvWriteActor._
	
	def receive = {
		case Write(traml, pcGroups, analysisTime) =>
			writeEsv(traml, pcGroups, analysisTime)
			sender ! FileWritten()
	}
	
	
	def writeEsv(
			traml:GhostTraML,
			pcGroups:Seq[Seq[DianaPeptideCandidate]], 
			analysisTime:Long
	) = {
		val before = System.currentTimeMillis
		
		var esv = new Esv
		esv.separator = '\t'
		esv.escape = '"'
		esv.nRows = Some(pcGroups.map(_.length).sum)
			
		esv.source = Some(params.name)
		esv.version = Some(params.version)
		
		esv.addParameter("upper bound"		, params.relRatioUpperBound)
		esv.addParameter("binSize"			, params.signalBinSize)
		esv.addParameter("p-value cutoff"	, params.pCutoff)
		esv.addParameter("swath file"		, params.swathFile)
		esv.addParameter("traml file"		, params.tramlFile)
		esv.addParameter("analysis time (ms)", analysisTime)
		
		/*if (AnalyzeCG.fragmentState.rtMapUsed) {
			esv.addParameter("rt map intercept", AnalyzeCG.fragmentState.rtIntersect)
			esv.addParameter("rt map slope", AnalyzeCG.fragmentState.rtSlope)
			esv.addParameter("rt map std", AnalyzeCG.fragmentState.rtStd)
		}
			*/
		
		esv.headers = Array("rawArea",
							"correctedArea",
							"isotopeArea",
						//	"rtProb",
							"fragmentRankAllRatioProb","fragmentRankPcsRatioProb","fragmentMarkovAllRatioProb","fragmentMarkovPcsRatioProb","fragmentCorrScore",
							"isotopeRankAllRatioProb","isotopeRankPcsRatioProb","isotopeMarkovAllRatioProb","isotopeMarkovPcsRatioProb","isotopeCorrScore",
							//"p-score","q-value",
							"rtStart","rtApex","rtEnd",
							"rtApexAssay",
							"q1","charge",
							"maxIntensity", "estimateApex",
							"alternatives",
							"peptideSequence","protein")
		esv.nColumns = Some(esv.headers.length)
		
		var bw = new BufferedWriter(new FileWriter(
					new File(params.outFile.toString + ".esv")
				))
		var ew = new EsvWriter(esv, bw)
		
		val uk = Locale.UK
		for (dpc <- pcGroups.flatten) {
			val (pep, prot) = 
				if (traml.peptides.contains(dpc.assay.pepCompRef)) {
					val pep = traml.peptides(dpc.assay.pepCompRef)
					(pep.sequence, pep.proteins.mkString)
				} else {
					(dpc.assay.pepCompRef, "-")
				}
			ew.write(List(
					"%10.1f".formatLocal(uk, dpc.rawArea), 
					"%10.1f".formatLocal(uk, dpc.correctedArea), 
					"%10.1f".formatLocal(uk, dpc.isotopeArea), 
					//"%.2e".formatLocal(uk, dpc.rtProb), 
					
					"%.2e".formatLocal(uk, dpc.fragmentRankAllRatioProb), 
					"%.2e".formatLocal(uk, dpc.fragmentRankPcsRatioProb), 
					"%.2e".formatLocal(uk, dpc.fragmentMarkovAllRatioProb), 
					"%.2e".formatLocal(uk, dpc.fragmentMarkovPcsRatioProb), 
					"%.4f".formatLocal(uk, dpc.fragmentCorrScore), 
					
					"%.2e".formatLocal(uk, dpc.isotopeRankAllRatioProb),
					"%.2e".formatLocal(uk, dpc.isotopeRankPcsRatioProb),
					"%.2e".formatLocal(uk, dpc.isotopeMarkovAllRatioProb),
					"%.2e".formatLocal(uk, dpc.isotopeMarkovPcsRatioProb),
					"%.4f".formatLocal(uk, dpc.isotopeCorrScore), 
					
					//"%.2e".formatLocal(uk, dpc.pscore), 
					//"%.2e".formatLocal(uk, dpc.qvalue), 
					"%.1f".formatLocal(uk, dpc.rtStart), 
					"%.1f".formatLocal(uk, dpc.rtApex), 
					"%.1f".formatLocal(uk, dpc.rtEnd), 
					"%5.1f".formatLocal(uk, dpc.assay.expectedRT),
					"%.2f".formatLocal(uk, dpc.assay.transitions.head.q1),
					dpc.assay.pepCompCharge,
					"%.1f".formatLocal(uk, dpc.maxIntensity),
					"%.1f".formatLocal(uk, dpc.maxEstimate),
					dpc.alternatives,
					pep, 
					prot
					))
		}
		
		ew.close
		params.esvWriteTime = System.currentTimeMillis - before
	}
}