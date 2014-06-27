package se.lth.immun.diana

import se.lth.immun.xml.XmlReader
import se.lth.immun.xml.XmlWriter
import se.lth.immun.app.CLIApplication
import se.lth.immun.app.CLIBar
import se.lth.immun.app.CommandlineArgumentException

import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.BufferedReader
import java.io.FileWriter
import java.io.BufferedWriter
import java.io.IOException
import java.util.Properties
import java.util.Date
import java.util.Locale

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashMap

import se.lth.immun.mzml._
import se.lth.immun.mzml.ghost._
import se.lth.immun.traml.ghost._

//import se.lth.immun.anubis.Ratio
//import se.lth.immun.anubis.AnubisInput
//import se.lth.immun.anubis.ReferenceFile
//import se.lth.immun.anubis.ReferencePrecursor
//import se.lth.immun.anubis.ReferenceTransition
//import se.lth.immun.anubis.AnubisParams
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

import se.lth.immun.esv._
import se.lth.immun.chem._
import se.lth.immun.unimod.UniMod

object DianaScorer extends CLIApplication {

	class Isotope(
			val q1:Double,
			val occurence:Double
	) {}
	
	
	// real
	var swathFile:File 			= _
	var tramlFile:File 			= _
	var traml:GhostTraML		= _
	//	var ref:ReferenceFile		= null
	var swathGmzML:GhostMzML	= _
	var outFile:File			= _
	var outEsv					= true
	var outXml					= false
	
	var noDecoys				= true

	// swath analysis
	var pScoreFunc:Carrier => Double = c => c.fragmentRankPcsRatioProb
	var nIsotopes				= 0
	
	var quiet 		= true
	var noBar		= false
	var cliBar		= new CLIBar(30, true)
	//var singleResult 	= true
	
	
	//def MISSING:Seq[PCGroup] = { var pg = new PCGroup; pg.pvalue = 1.0; Array(pg) }
	
	
	def main(args:Array[String]):Unit = {
		
		var properties = new Properties
    	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
    	
    	arg("SWATH_CHROM_MZML", s => {
			swathFile = new File(s)
			if (outFile == null) {
				val comp = swathFile.toString.toLowerCase
				if (comp.endsWith(".chrom.mzml")) 
					outFile = new File(swathFile.toString.dropRight(11))
				else if (comp.endsWith(".mzml"))
					outFile = new File(swathFile.toString.dropRight(5))
				else
					outFile = new File(swathFile.toString)
			}
		})
		
    	arg("TRAML", s => {
			tramlFile = new File(s)
			val ext = s.split('.').last.toLowerCase
			if (ext == "traml")
				traml = GhostTraML.fromFile(new XmlReader(
												new BufferedReader(new FileReader(tramlFile))
											))
			/*	ref = ReferenceFile.fromTraML(new XmlReader(
												new BufferedReader(new FileReader(tramlFile))
											))
			else if (ext == "ref" || ext == "xml")
				ref = ReferenceFile.fromFile(new XmlReader(
												new BufferedReader(new FileReader(tramlFile))
											))
			*/
			else {
				val s = "Unkown traml extention '%s'. Need .traml".format(ext) //, .ref och .xml".format(ext)
				println(s)
				throw new IllegalArgumentException(s)
			}
		})
		
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
		
		/*opt("no-bar", 
				"don't print process bar", 
				s => noBar = true)
		
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
		*/
		opt("verbose", 
				"(default: not verbose)", 
				s => {
					quiet = false
				})
		
				/*
		opt("single-result", 
				"report more than 1 peak per assay if good enough peaks exist (default: not used)", 
				s => {
					singleResult = false
				})
		*/
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
		
				/*
		
		opt("isotopes", 
				"The number highest precursor isotopes to analyze (default 0)", 
				s => nIsotopes = s.toInt, "X")
		*/
		
		val before 		= System.currentTimeMillis
    	val name 		= properties.getProperty("pom.name")
    	val version 	= properties.getProperty("pom.version")
    	
	/*
	val f = new java.text.DecimalFormatSymbols(Locale.getDefault())
	println("java dec sep: "+f.getDecimalSeparator)
	println("scala float format:   %.3f %.3f".format(0.029, 23000.0))
	println("scala science format: %.3e %.3e".format(0.029, 23000.0))
	println("scala float format UK:   %.3f %.3f".formatLocal(Locale.UK, 0.029, 23000.0))
	println("scala science format UK: %.3e %.3e".formatLocal(Locale.UK, 0.029, 23000.0))
	*/
		try {
			parseArgs(name + " "+version, args)
		} catch {
			case cae:CommandlineArgumentException => {
				cae.printStackTrace
				return
			}
		}

					
							
		println(name + " "+version)
    	println("  swath mzML file: " + swathFile)
		println("       traML file: " + tramlFile)
		println("         p cutoff: " + AnalyzeCG.pValueCutoff)
    	println("     num isotopes: " + nIsotopes)
		//if (noDecoys)
		//	println("                => no decoy information, will skip FDR calculations")
    	println()
		
		var qvalues:Seq[(Double, Double)] = Nil
    	var pcgroups = analyzeFile(swathFile, traml, false)
    					.map(_.sortBy(_.pscore))
    					.sortBy(_.head.pscore)
    		/*
    	if (singleResult)
    		pcgroups = pcgroups.map(_.headOption.toSeq)
		
		*/
		val analysisTime = System.currentTimeMillis - before
		/*if (outXml) {
			var res = new ResultFile(Array(),
					Array(new ResultMzMLFile(swathFile, swathGmzML.runTime)))
			
			res.referenceFile = tramlFile
			res.workingDir = new File(".")
			
			//println("Q1 pval qval area rt")
			res.precursors = pcgroups.flatten.map(dpc => {
				val assay = dpc.assay
				
				//println("%6.1f %10.4e %10.4e %10.0f %4.1f".format(rp.mz, spc.pvalue, spc.qvalue, spc.area, spc.rtApex))
				new ResultPrecursor(
					assay.transitions.head.q1,
					rp.proteinName,
					rp.peptideSequence,
					Array(
						new ResultReplicate(
							spc.correctedArea,
							0,
							if (!spc.missing)
								rp.measuredTransitions.zip(spc.correctedAreas).map(t => 
									new ResultTransition(
										t._1.mz,
										t._1.ion,
										t._2
									)
								).toArray
							else 
								rp.measuredTransitions.map(t => 
									new ResultTransition(
										t.mz,
										t.ion,
										0.0
									)
								).toArray,
							new ResultRetentionTime(
								spc.rtStart, spc.rtApex, spc.rtEnd
							),
							new ResultQuality(if (spc.missing) Double.NaN else if (noDecoys) (1-spc.pscore) else spc.qvalue)
						)
					)
				)
			}).toArray
			
			res.write(new XmlWriter(new BufferedWriter(new FileWriter(
						new File(outFile.toString + ".xml")
					))))
		}
		*/
		if (outEsv) {
			var esv = new Esv
			esv.separator = '\t'
			esv.escape = '"'
			esv.nRows = Some(pcgroups.map(_.length).sum)
				
			esv.source = Some(name)
			esv.version = Some(version)
			
			esv.addParameter("upper bound", AnalyzeCG.fragmentState.upperBound)
			esv.addParameter("binSize", AnalyzeCG.fragmentState.binSize)
			esv.addParameter("p-value cutoff", AnalyzeCG.pValueCutoff)
			esv.addParameter("used decoys", !noDecoys)
			//esv.addParameter("single result", singleResult)
			esv.addParameter("swath file", swathFile)
			esv.addParameter("traml file", tramlFile)
			esv.addParameter("analysis time (ms)", analysisTime)
			/*
			if (decoySwathFile != null) {
				esv.addParameter("decoy swath file", decoySwathFile)
				esv.addParameter("decoy traml file", decoyTramlFile)
			}
			if (decoyEsvFile != null)
				esv.addParameter("decoy esv file", decoyEsvFile)
			 */
			if (AnalyzeCG.fragmentState.rtMapUsed) {
				esv.addParameter("rt map intercept", AnalyzeCG.fragmentState.rtIntersect)
				esv.addParameter("rt map slope", AnalyzeCG.fragmentState.rtSlope)
				esv.addParameter("rt map std", AnalyzeCG.fragmentState.rtStd)
			}
				
			
			esv.headers = Array("rawArea",
								"correctedArea",
								"isotopeArea",
								"rtProb",
								"fragmentRankAllRatioProb","fragmentRankPcsRatioProb","fragmentMarkovAllRatioProb","fragmentMarkovPcsRatioProb","fragmentCorrScore",
								"isotopeRankAllRatioProb","isotopeRankPcsRatioProb","isotopeMarkovAllRatioProb","isotopeMarkovPcsRatioProb","isotopeCorrScore",
								"p-score","q-value",
								"rtStart","rtApex","rtEnd",
								"rtApexAssay",
								"q1","charge",
								"maxIntensity", "estimateApex",
								"alternatives",
								"protein","peptideSequence")
			esv.nColumns = Some(esv.headers.length)
			
			var bw = new BufferedWriter(new FileWriter(
						new File(outFile.toString + ".esv")
					))
			var ew = new EsvWriter(esv, bw)
			
			val uk = Locale.UK
			for (dpc <- pcgroups.flatten) {
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
						"%.2e".formatLocal(uk, dpc.rtProb), 
						
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
						
						"%.2e".formatLocal(uk, dpc.pscore), 
						"%.2e".formatLocal(uk, dpc.qvalue), 
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
		}
		val after = System.currentTimeMillis
		println("  time taken: "+(after-before)+"ms")
		
		return
	}
	
	
	
	def analyzeFile(
			mzmlFile:File,
			traml:GhostTraML,
			decoy:Boolean
	):Seq[Seq[DianaPeptideCandidate]] = {
		val (reader, binFile, binFileChannel) = GhostMzML.getReaders(mzmlFile)
		var gmzML = GhostMzML.fromFile(reader, true, binFileChannel)
		swathGmzML = gmzML
		
		val assays = getAssays(traml)
		
		var aCount 	= 0
		val nAssays 	= assays.length
		cliBar.reset
		if (quiet) {
			println(cliBar.reference)
			print(cliBar.update(0))
		}
		
		assays.map(assay => {
			if (quiet) {
				aCount += 1
				print(cliBar.update(aCount, nAssays))
			}
			val transChroms = assay.transitions.map(t => gmzML.chromatograms.find(_.id == t.id))
			val targChroms = assay.targets.map(t => gmzML.chromatograms.find(_.id == t.id))
			
			if (transChroms.count(_.isDefined) < 2) {
				val errMsg = "ERROR during analysis of '%s %.3f' in file '%s'".format(assay.pepCompRef, assay.pepCompMz, mzmlFile.toString)
        		CLIApplication.log.write(errMsg)
            	if (!quiet) 
            		println(errMsg)
            	List() //List(new DianaPeptideCandidate(pc))
			} else {
				
				val dIn = DianaInput.fromTransitions(transChroms.map(_.get.toXChromatogram), assay)
				
				val iDIn =
					if (targChroms.nonEmpty)
						Some(DianaInput.fromTargets(
								targChroms.map(_.get.toXChromatogram), 
								assay, dIn.cg.chromatograms.head.times
							))
					else None
					
				val carriers = AnalyzeCG(dIn, iDIn, assay)
				if (carriers.isEmpty) List(new DianaPeptideCandidate(assay))
		        else carriers.map(c => DianaPeptideCandidate(assay, dIn.cg, iDIn.map(_.cg), c, pScoreFunc, carriers.length))
		       /* 
				try {
		        	xmzML.grouper.extractGroup(pc.mz) match {
		            	case Some(cg) => {
		            		val aIn = AnubisInput.of(cg, pc, new AnubisParams)
		            		
		            		var iAIn =
		            			if (nIsotopes > 1)
		            				getIsotopeInput(xmzML, pc, aIn.cg.chromatograms.head.times)
		            			else None
		            		
		            		val carriers = AnalyzeCG(aIn, iAIn, pc)
		            		if (carriers.isEmpty) List(new DianaPeptideCandidate(pc))
		            		else carriers.map(c => DianaPeptideCandidate(pc, aIn.cg, iAIn.map(_.cg), c, pScoreFunc, carriers.length))
		            	}
		            	case None => {
		            		CLIApplication.log.write("%s %.3f not found in '%s'".format(pc.peptideSequence, pc.mz, mzmlFile.toString))
		            		List(new DianaPeptideCandidate(pc))
		            	}
		            }
		        }
		        catch {
		        	case e:Exception => {
		        		val errMsg = "ERROR during analysis of '%s %.3f' in file '%s'".format(pc.peptideSequence, pc.mz, mzmlFile.toString)
		        		CLIApplication.log.write(errMsg)
		            	CLIApplication.log.write(e)
		            	if (!quiet) {
		            		println(errMsg)
		            		e.printStackTrace
		            	}
		        		List(new DianaPeptideCandidate(pc))
		        	}
		        }*/
			}
		})
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
							CLIApplication.log.write(errMsg)
							if (!quiet) 
								println(errMsg)
							0
					}
				} else {
					val errMsg = "WARN: no charge given for compound '%s' with mz %.3f".format(t.compoundRef, t.q1)
					CLIApplication.log.write(errMsg)
					if (!quiet) 
						println(errMsg)
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
							CLIApplication.log.write(errMsg)
							if (!quiet) 
								println(errMsg)
							0
					}
				} else {
					val errMsg = "WARN: no charge given for compound '%s' with mz %.3f".format(t.compoundRef, t.q1)
					CLIApplication.log.write(errMsg)
					if (!quiet) 
						println(errMsg)
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
