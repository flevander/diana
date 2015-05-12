package se.lth.immun.diana

import se.lth.immun.xml.XmlReader
import se.lth.immun.xml.XmlWriter
import se.jt.CLIApp

import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.BufferedReader
import java.io.IOException
import java.util.Properties
import java.util.Date

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashMap

import se.lth.immun.mzml._
import se.lth.immun.mzml.ghost._
import se.lth.immun.traml.ghost._

import se.lth.immun.anubis.ResultFile
import se.lth.immun.anubis.ResultParameter
import se.lth.immun.anubis.ResultMzMLFile
import se.lth.immun.anubis.ResultPrecursor
import se.lth.immun.anubis.ResultReplicate
import se.lth.immun.anubis.ResultTransition
import se.lth.immun.anubis.ResultRetentionTime
import se.lth.immun.anubis.ResultQuality

import DianaPeakCandidate._

import se.lth.immun.math.Ratios

import se.lth.immun.chem._
import se.lth.immun.unimod.UniMod

import akka.actor.Actor
import akka.actor.ActorRef
import akka.actor.ActorSystem
import akka.actor.Props
import akka.actor.Inbox

object DianaScorer extends CLIApp {

	class Isotope(
			val q1:Double,
			val occurence:Double
	) {}

	
	def main(args:Array[String]):Unit = {
		
		var properties = new Properties
    	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
    	val name 		= properties.getProperty("pom.artifactId")
    	val version 	= properties.getProperty("pom.version")
    	
    	val params = new DianaScorerParams(name, version)
    	
		params.startTime 	= System.currentTimeMillis
    	
    	failOnError(parseArgs(name, version, args, params, List("swathChromMzML", "traml"), None))
    	
		println(name + " "+version)
    	println("  swath mzML file: " + params.swathChromMzML.value)
		println("       traML file: " + params.traml.value)
		println()
		println(" reading traml...")
		
		
    	fixDefaultOutputPath(params)
    		
    	val (traml, errs) = readTraml(params)
		failOnError(errs)
    	
		import AnalysisActor._
	
    	val system = ActorSystem("diana-system")
		val top = system.actorOf(Props(new TopActor(params)), "top-actor")
    	top ! AnalyzeFile(traml, params)
    	
    	system.awaitTermination
    	
		/*
		opt("output", 
				"file where result should be saved (default: input chrom .xml/.esv)", 
				s => {
					outFile = new File(s)
				},
				"X")
		
		opt("p-cutoff", 
				"p-value cutoff for peak candidates (default: 0.99)", 
				s => {
					AnalyzeCG.pValueCutoff = s.toDouble
				},
				"X")
		
		opt("concurrency", 
				"the maximal amount of assays to compute in parallel (default: 1)", 
				s => {
					concurrency = s.toInt
				},
				"X")
		
		opt("edges", 
				"edge find strategy to use: strict, L1 och L2 (default: strict)", 
				s => {
					s match {
						case "L1" => findEdges = looseFindEdges(1) _
						case "L2" => findEdges = looseFindEdges(2) _
						case _ => {}
					}
				},
				"strict|L1|L2")
		
		opt("verbose", 
				"(default: not verbose)", 
				s => {
					quiet = false
				})
				
		opt("cpu-profiling", 
				"(default: no profiling)", 
				s => {
					profiling = true
				})
		
		opt("no-esv", 
				"flag if esv output it not wanted", 
				s => {
					outEsv = false
				})
		
		opt("xml", 
				"flag if xml output is wanted", 
				s => {
					outXml = true
				})
		
    	
		try {
			parseArgs(name + " "+version, args)
		} catch {
			case cae:CommandlineArgumentException => {
				cae.printStackTrace
				return
			}
		}
		 */
					
    	
	}
	
	
	
	
	def fixDefaultOutputPath(p:DianaScorerParams) = {
		p.swathFile = new File(p.swathChromMzML)
		if (p.outFile == null) {
			val comp = p.swathFile.toString.toLowerCase
			if (comp.endsWith(".chrom.mzml")) 
				p.outFile = new File(p.swathFile.toString.dropRight(11))
			else if (comp.endsWith(".mzml"))
				p.outFile = new File(p.swathFile.toString.dropRight(5))
			else
				p.outFile = new File(p.swathFile.toString)
		}
	}
	
	
	
	
	def readTraml(p:DianaScorerParams):(GhostTraML, List[String]) = {
		p.tramlFile = new File(p.traml)
		val before = System.currentTimeMillis
		val ext = p.traml.value.split('.').last.toLowerCase
		val traml =
			if (ext == "traml")
				GhostTraML.fromFile(new XmlReader(
											new BufferedReader(new FileReader(p.tramlFile))
										))
			else 
				return (null, List("Unkown traml extention '%s'. Need .traml".format(ext)))
		
		p.tramlReadTime = System.currentTimeMillis - before
		(traml, Nil)
	}
	
	
	
	
	
	
	
	/*
	def getIsotopeInput(xmzML:XMzML, pc:ReferencePrecursor, times:Seq[Double]):Option[AnubisInput] = {
		val seq 	= pc.peptideSequence
		val p 		= UniMod.parseUniModSequence(seq)
		val q1		= pc.mz
		val q1z 	= math.round(p.monoisotopicMass() / q1).toDouble
		val id = p.getIsotopeDistribution()
		val isotopes = new ArrayBuffer[Isotope]
		for (i <- id.intensities.sorted.takeRight(nIsotopes)) {
			val ii = id.intensities.indexOf(i)
			isotopes += new Isotope(q1 + ii / q1z, i)
		}
		
		val isotopeChroms = isotopes.map(i => { 
			var cgo = xmzML.grouper.extractGroup(i.q1, 0.01)
			cgo match {
				case Some(cg) => cg.chromatograms.find(_.q3 == 0.0).map(c => (i, c))
				case None => None
			}
		}).filter(_.isDefined).map(_.get)
		
		if (isotopeChroms.isEmpty || isotopeChroms.forall(_._2.intensities.max == 0)) 
			return None
		
		val cg = new XChromatogramGroup
		for (ic <- isotopeChroms)
			cg.chromatograms += ic._2
		
		val r = new ArrayBuffer[Ratio]
		Ratios.iterate(isotopeChroms.length, (ri, ui, di) => {
			r += new Ratio(ui, di, isotopeChroms(ui)._1.occurence / isotopeChroms(di)._1.occurence, -1.0)
		})
		
		return Some(new AnubisInput(
				seq, 
				q1, 
				cg.resample(times), 
				isotopeChroms.map(ic => {
						val mz = ic._1.q1
						AnubisInput.Fragment(mz, "isotope %.2f".formatLocal(Locale.UK, mz))
					}), 
				r.toArray, 
				new AnubisParams
			))
	}
	
	*/
	
	
	
	
	
	
	/*
	def calculateQvalues(
			groups:Seq[Seq[DianaPeptideCandidate]], 
			decoys:Seq[Double]
	):Seq[Seq[DianaPeptideCandidate]] = {
		var ig = 0
		var id = 0
		var gdratio = groups.length.toDouble / decoys.length
		var maxq = 1.0 / (groups.length * decoys.length)
		while (ig < groups.length && id < decoys.length) {
			if (groups(ig).head.pscore < decoys(id)) {
				maxq = math.max(
						maxq, 
						(id.toDouble / (ig+1)) * gdratio
						)
				groups(ig).head.qvalue = maxq
				ig += 1
			} else {
				id += 1
			}
		}
		for (i <- ig until groups.length)
			groups(i).head.qvalue = 1.0
		
		val pqscale = groups.map(s => (s.head.pscore, s.head.qvalue))
		for (g <- groups)
			if (g.length > 1)
				for (spc <- g.tail)
					pqscale.find(t => t._1 > spc.pscore) match {
						case Some(t) => 	spc.qvalue = t._2
						case None => 		spc.qvalue = 1.0
					}
		
		return groups
	}
	*/
}
