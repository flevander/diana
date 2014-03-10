package se.lth.immun.diana

import se.lth.immun.xml.XmlReader
import se.lth.immun.xml.XmlWriter
import se.lth.immun.app.CLIApplication
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

import se.lth.immun.mzml._
import se.lth.immun.mzml.ghost._
import se.lth.immun.traml.ghost._

import se.lth.immun.anubis.Ratio
import se.lth.immun.anubis.AnubisInput
import se.lth.immun.anubis.ReferenceFile
import se.lth.immun.anubis.ReferencePrecursor
import se.lth.immun.anubis.ReferenceTransition
import se.lth.immun.anubis.AnubisParams
import se.lth.immun.anubis.ResultFile
import se.lth.immun.anubis.ResultParameter
import se.lth.immun.anubis.ResultMzMLFile
import se.lth.immun.anubis.ResultPrecursor
import se.lth.immun.anubis.ResultReplicate
import se.lth.immun.anubis.ResultTransition
import se.lth.immun.anubis.ResultRetentionTime
import se.lth.immun.anubis.ResultQuality

import SwathPeakCandidate._

import se.lth.immun.math.Ratios

import se.lth.immun.esv._
import se.lth.immun.chem._

object DianaScorer extends CLIApplication {

	class Isotope(
			val q1:Double,
			val occurence:Double
	) {}
	
	
	// real
	var swathFile:File 			= null
	var tramlFile:File 			= null
	var ref:ReferenceFile		= null
	var swathXmzML:XMzML		= null
	var outFile:File			= null
	var outEsv					= true
	var outXml					= true
	// decoy
	/*
	var decoySwathFile:File 	= null
	var decoyTramlFile:File 	= null
	var decoyRef:ReferenceFile	= null
	var decoyXmzML:XMzML		= null
	var decoyEsvFile:File		= null
	var decoyEsv:EsvReader		= null
	var decoyPValues:Seq[Double] = Nil
		*/
	var noDecoys				= true

	// swath analysis
	var pScoreFunc:Carrier => Double = c => c.fragmentRankPcsRatioProb
	var nIsotopes				= 0
	
	var quiet 			= true
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
				ref = ReferenceFile.fromTraML(new XmlReader(
												new BufferedReader(new FileReader(tramlFile))
											))
			else if (ext == "ref" || ext == "xml")
				ref = ReferenceFile.fromFile(new XmlReader(
												new BufferedReader(new FileReader(tramlFile))
											))
			else {
				val s = "Unkown traml extention '%s'. Need .traml, .ref och .xml".format(ext)
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
		/*
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
				/*
		opt("rt-prob", 
				"file w. expected empirical linear rt model vs reference (default: rtProb = 0)", 
				s => {
					var r = new BufferedReader(new FileReader(new File(s)))
					r.readLine
					val t = r.readLine.split(" ")
					AnalyzeCG.fragmentState.rtIntersect = t(0).toDouble
					AnalyzeCG.fragmentState.rtSlope 	= t(1).toDouble
					AnalyzeCG.fragmentState.rtStd 		= t(2).toDouble
				},
				"RT_FILE")
		
		opt("p-score", 
				"file w. linear model for computing p-score (default: p-score = nullRatioProb)", 
				s => {
					var r = new BufferedReader(new FileReader(new File(s)))
					r.readLine
					val t = r.readLine.split(" ")
					val intercept 	= t(0).toDouble
					val rtK 		= t(1).toDouble
					val ratioK 		= t(2).toDouble
					val fragK 		= t(3).toDouble
					val corrK 		= t(4).toDouble
					pScoreFunc = g => 
						intercept + g.rtProb*rtK + g.fragmentRatioProb*ratioK + g.pEstimates*fragK + g.fragmentCorrScore*corrK
				},
				"P_SCORE_FILE")
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
		
		opt("no-xml", 
				"flag if xml output it not wanted", 
				s => {
					outXml = false
				})
		
				/*
		opt("decoy-traml", 
				"Traml with decoy transitions, also need --decoy-mzml with the corresponding data", 
				s => {
					decoyTramlFile = new File(s)
					decoyRef = ReferenceFile.fromTraML(new XmlReader(
											new BufferedReader(new FileReader(decoyTramlFile))
										))
				},
				"TRAML"
		)
		
		opt("decoy-mzml", 
				"Mzml with decoy data, also need --decoy-traml with the corresponding transitions", 
				s => {
					decoySwathFile = new File(s)
				},
				"MZML"
		)
		
		opt("decoy-esv", 
				"esv with decoy data, here the p-values are precalculated", 
				s => {
					decoyEsvFile = new File(s)
					decoyEsv = new EsvReader(new BufferedReader(new FileReader(decoyEsvFile)))
				},
				"ESV"
		)
		*/
		opt("isotopes", 
				"The number highest precursor isotopes to analyze (default 0)", 
				s => nIsotopes = s.toInt, "X")
		
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

		
		/*
		if (decoyEsv != null) {
			if (decoyTramlFile != null || decoySwathFile != null)
				return println("Argument error: both decoy-esv and one or both of decoy-mzml and decoy-traml set. "+
							"Use only decoy-esv or both of decoy-mzml and decoy-traml. Exiting.")
			noDecoys = false
		} else {
			if (decoyTramlFile != null && decoySwathFile == null)
				return println("Argument error: decoy-traml set but not decoy-mzml. "+
								"Need both for decoy handling. Exiting.")
			if (decoyTramlFile == null && decoySwathFile != null)
				return println("Argument error: decoy-mzml set but not decoy-traml. "+
								"Need both for decoy handling. Exiting.")
			if (decoyTramlFile != null && decoySwathFile != null)
				noDecoys = false
		}
			*/				
							
		println(name + " "+version)
    	println("  swath mzML file: " + swathFile)
		println("       traML file: " + tramlFile)
		println("         p cutoff: " + AnalyzeCG.pValueCutoff)
    	println("     num isotopes: " + nIsotopes)
		//if (noDecoys)
		//	println("                => no decoy information, will skip FDR calculations")
    	println()
    	
    	
		var qvalues:Seq[(Double, Double)] = Nil
    	var pcgroups = analyzeFile(swathFile, ref, false)
    					.map(_.sortBy(_.pscore))
    					.sortBy(_.head.pscore)
    		/*
    	if (singleResult)
    		pcgroups = pcgroups.map(_.headOption.toSeq)
		
		if (!noDecoys) {
    		if (decoyEsv != null) {
    			var a = new ArrayBuffer[Double]
    			while (!decoyEsv.EOF) {
    				a += decoyEsv.getValue("p-score").toDouble
    				decoyEsv.readLine
    			}
    			decoyEsv.close
    			decoyPValues = a
    		} else {
    			decoyPValues = analyzeFile(decoySwathFile, decoyRef, true)
    				.map(_.map(_.pscore).min).sorted
    		}
    		
    		if (!quiet) {
    			println("     REAL: "+pcgroups.map(c => "%.4f".format(c.head.pscore)))
    			println("    DECOY: "+decoyPValues.map(pscore => "%.4f".format(pscore)))
    		}
    		qvalues = calculateQvalues(pcgroups, decoyPValues)
    					.map(x => (x.head.pscore, x.head.qvalue))
    		
    		if (!quiet)
    			println(" Q-VALUES: "+qvalues.map(t => 
    								"%.4f".format(t._2)
    							).mkString(" "))
    	}
		*/
		val analysisTime = System.currentTimeMillis - before
		if (outXml) {
			var res = new ResultFile(Array(),
					Array(new ResultMzMLFile(swathFile, swathXmzML.runTime)))
			/*
			var res = new ResultFile(
					if (noDecoys) Array()
					else if (decoyEsvFile != null)
						Array(new ResultParameter("decoy esv file", decoyEsvFile.toString))
					else Array(
						new ResultParameter("decoy mzML file", decoySwathFile.toString),
						new ResultParameter("decoy traML file", decoyTramlFile.toString)
					),
					Array(new ResultMzMLFile(swathFile, swathXmzML.runTime))
				)
	    	*/
			res.referenceFile = tramlFile
			res.workingDir = new File(".")
			
			//println("Q1 pval qval area rt")
			res.precursors = pcgroups.flatten.map(spc => {
				val rp = spc.reference
				//println("%6.1f %10.4e %10.4e %10.0f %4.1f".format(rp.mz, spc.pvalue, spc.qvalue, spc.area, spc.rtApex))
				new ResultPrecursor(
					rp.mz,
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
			for (spc <- pcgroups.flatten) 
				ew.write(List(
						"%10.1f".formatLocal(uk, spc.rawArea), 
						"%10.1f".formatLocal(uk, spc.correctedArea), 
						"%.2e".formatLocal(uk, spc.rtProb), 
						
						"%.2e".formatLocal(uk, spc.fragmentRankAllRatioProb), 
						"%.2e".formatLocal(uk, spc.fragmentRankPcsRatioProb), 
						"%.2e".formatLocal(uk, spc.fragmentMarkovAllRatioProb), 
						"%.2e".formatLocal(uk, spc.fragmentMarkovPcsRatioProb), 
						"%.4f".formatLocal(uk, spc.fragmentCorrScore), 
						
						"%.2e".formatLocal(uk, spc.isotopeRankAllRatioProb),
						"%.2e".formatLocal(uk, spc.isotopeRankPcsRatioProb),
						"%.2e".formatLocal(uk, spc.isotopeMarkovAllRatioProb),
						"%.2e".formatLocal(uk, spc.isotopeMarkovPcsRatioProb),
						"%.4f".formatLocal(uk, spc.isotopeCorrScore), 
						
						"%.2e".formatLocal(uk, spc.pscore), 
						"%.2e".formatLocal(uk, spc.qvalue), 
						"%.1f".formatLocal(uk, spc.rtStart), 
						"%.1f".formatLocal(uk, spc.rtApex), 
						"%.1f".formatLocal(uk, spc.rtEnd), 
						"%5.1f".formatLocal(uk, spc.reference.retentionTime.peak),
						"%.2f".formatLocal(uk, spc.reference.mz),
						spc.reference.charge,
						"%.1f".formatLocal(uk, spc.maxIntensity),
						"%.1f".formatLocal(uk, spc.maxEstimate),
						spc.alternatives,
						spc.reference.proteinName, 
						spc.reference.peptideSequence
						))
			
			ew.close
		}
		val after = System.currentTimeMillis
		println("  time taken: "+(after-before)+"ms")
		
		return
	}
	
	
	
	def analyzeFile(
			mzmlFile:File,
			ref:ReferenceFile,
			decoy:Boolean
	):Seq[Seq[SwathPeptideCandidate]] = {
		val (reader, binFile, binFileChannel) = XMzML.getReaders(mzmlFile)
		var xmzML = XMzML.fromFile(reader, true, binFileChannel)
		swathXmzML = xmzML
		/*if (decoy) {
			decoyXmzML = xmzML
		} else {
			swathXmzML = xmzML
		}
		*/
		ref.precursors.map(pc => {
			try {
	        	xmzML.grouper.extractGroup(pc.mz) match {
	            	case Some(cg) => {
	            		var iAIn:AnubisInput = null
	            		if (nIsotopes > 1) {
	            			iAIn = getIsotopeInput(xmzML, pc)
	            		}
	            		
	            		val aIn = AnubisInput.of(cg, pc, new AnubisParams)
	            		val carriers = AnalyzeCG(aIn, iAIn, pc)
	            		if (carriers.isEmpty) List(new SwathPeptideCandidate(pc))
	            		else carriers.map(c => SwathPeptideCandidate(pc, aIn.cg, c, pScoreFunc, carriers.length))
	            	}
	            	case None => {
	            		//TODO: error handling
	            		List(new SwathPeptideCandidate(pc))
	            	}
	            }
	        }
	        catch {
	        	case e:Exception => {
	            	//TODO: error handling
	        		e.printStackTrace
	        		List(new SwathPeptideCandidate(pc))
	        	}
	        }
		})
	}
	
	
	
	def getIsotopeInput(xmzML:XMzML, pc:ReferencePrecursor):AnubisInput = {
		val seq 	= pc.peptideSequence
		val p 		= new Peptide(seq.map(c => StandardAminoAcid.fromChar(c)).filter(_ != null).toArray)
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
		
		if (isotopeChroms.isEmpty) return null
		val cg = new XChromatogramGroup
		for (ic <- isotopeChroms)
			cg.chromatograms += ic._2
		
		val r = new ArrayBuffer[Ratio]
		Ratios.iterate(isotopeChroms.length, (ri, ui, di) => {
			r += new Ratio(ui, di, isotopeChroms(ui)._1.occurence / isotopeChroms(di)._1.occurence, -1.0)
		})
		
		return new AnubisInput(seq, q1, cg, isotopeChroms.map(ic => {
						val mz = ic._1.q1
						AnubisInput.Fragment(mz, "isotope %.2f".formatLocal(Locale.UK, mz))
					}), 
					r.toArray, new AnubisParams)
	}
	
	
	
	
	
	
	
	
	/*
	def calculateQvalues(
			groups:Seq[Seq[SwathPeptideCandidate]], 
			decoys:Seq[Double]
	):Seq[Seq[SwathPeptideCandidate]] = {
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
