package se.lth.immun.diana

import se.lth.immun.xml.XmlReader
import se.lth.immun.xml.XmlWriter
import se.jt.CLIApp
import se.jt.CLIArgumentException

import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.BufferedReader
import java.io.FileWriter
import java.io.BufferedWriter
import java.io.IOException
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream
import java.util.Calendar
import java.util.Properties
import java.text.SimpleDateFormat

import se.lth.immun.mzml._
import se.lth.immun.mzml.ghost._
import se.lth.immun.traml.ghost._
import se.lth.immun.chem._
import se.lth.immun.unimod.UniMod

import ms.numpress.MSNumpress

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.ArrayBuilder

import org.apache.commons.math3.distribution.NormalDistribution

object DianaExtractor extends CLIApp {

	import MzML._
	
	val SWATHS_IN_FILE = 32
	object CHROM_FILE_CONTENT extends FileContent {
		cvParams += {
			val x = new CvParam
			x.cvRef = "MS"
			x.accession = "MS:1001473"
			x.name = "selected reaction monitoring chromatogram"
			x
		}
	}
	
	object SOFTWARE extends Software {
		id = "SWATHExtractor"
		version = "0.2.3-SNAPSHOT"
	}
	
	object DATA_PROCESSING extends DataProcessing {
		
		var cv = new se.lth.immun.mzml.CvParam
		cv.cvRef 		= "MS"
		cv.accession 	= "MS:1000035"
		cv.name 		= "chromatogram extraction"
		
		var pm = new ProcessingMethod
		pm.order = 1
		pm.softwareRef = SOFTWARE.id
		pm.cvParams += cv
		
		id = "SWATHExtraction"
		processingMethods += pm
	}
	
	class Isotope(
			val seq:String,
			val q1:Double,
			val occurence:Double
	) {}
	
	
	
	
	val MODE_BEST = 0
	val MODE_UNIFORM = 1
	val MODE_NORMAL = 2
	val MODE_SQUARE = 3
	val MODE_TRI = 4
	
	
	var tramls:Seq[GhostTraML]	= Nil
	var xr:XmlReader 			= _
	var numSpec					= 0
	var specPerSwath			= 0
	var getQ1Window:(GhostSpectrum => TransitionPartitioning.Q1Window) = gilletQ1Window _
	var swathWidth 				= 25.0
	var mode					= MODE_UNIFORM
	val nd 					= new NormalDistribution

	var transitionSets		= new ArrayBuffer[TransitionPartitioning]
	var isotopeSets			= new ArrayBuffer[Seq[ChromExtract[Isotope]]]
	//var isotopeChromSets	= new ArrayBuffer[Array[Array[Double]]]
	var times 				= new Array[Array[Double]](SWATHS_IN_FILE)
	var ms1Times:ChromBuilder = _
	//var specCounters		= new Array[Int](SWATHS_IN_FILE)
	
	var ms1SpecCounter		= -1
	var capacitiesFixed		= false
	var capFixPercent		= 0.1
	
	
	var outFile:File 	= null
	var out:XmlWriter 	= null
	var currTramlOut	= -1
	
	var barCount = 0
	var specCount = 0
	
	def isIMzML(f:File) 	= f.getName.toLowerCase.endsWith(".imzml")
	def toIBD(f:File)		= new File(f.getAbsolutePath().dropRight(5) + "ibd")
	
	var properties = new Properties
	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
	val name 		= properties.getProperty("pom.name")
	val version 	= properties.getProperty("pom.version")
	val params = new DianaExtractorParams
	
	def main(args:Array[String]):Unit = {
				
		val t0 = System.currentTimeMillis
    	
    	failOnError(parseArgs(name, version, args, params, List("mzML"), Some("tramls")))
    	
		
		
		if (!params.mzMLFile.exists)
			failOnError(List("Input mzML file '%s' does not exists!".format(params.mzML.value)))
			
		tramls = params.tramlFiles.map(f => 
				GhostTraML.fromFile(
					new XmlReader(
						new BufferedReader(new FileReader(f))
					)
				)
			)
		
		println(name + " "+version)
    	println("dia/swath mzML file: " + params.mzML.value)
    	println("       num isotopes: " + params.nIsotopes.value)
    	println("     sub sample ms1: " + params.subSampleMs1.value)
    	println("              force: " + params.force.value)
    	println("      ms resolution: " + params.mzThreshold)
    	println("    extraction mode: " + params.mode.value)
    	println("   dia-windows mode: " + params.diaWindows.value)
		println("        traML files:")
    		params.tramlFiles.foreach(f => println("    "+f))
    	println()
    	
		xr = params.mzMLReader
		xr.force = params.force
		handleSwathFile(xr, params.mzMLFile)
		 
		val t1 = System.currentTimeMillis
		val nSpectraPerPartition = transitionSets.head.parts.values.map(_.times.nocopyResult.length)
		val runt = Runtime.getRuntime
		println("        heap size: "+(runt.totalMemory / 1000000) + " / "+(runt.maxMemory / 1000000) + " Mb")
		println("    n ms1 spectra: "+ms1SpecCounter)
		println("    n dia windows: "+transitionSets.head.parts.size)
		println("spectra per swath: "+nSpectraPerPartition.min + " - "+nSpectraPerPartition.max)
		println("       time taken: "+niceTiming(t1-t0))
		println("     spectra read: "+transitionSets.head.parts.values.map(_.times.nocopyResult.length).sum)
	}
	
	
	def parseDiaWindows(s:String) = 
		if (Array("gillet", "snoop").contains(s)) {
			if (s == "gillet") 		getQ1Window = gilletQ1Window _
			else if (s == "snoop") 	getQ1Window = snoopQ1Window _
		} else 
			throw new Exception("Unknown dia window mode '"+s+"'")
		
	def parseMode(s:String) =
		s match {
			case "best" 	=> MODE_BEST
			case "uniform" 	=> MODE_UNIFORM
			case "normal" 	=> MODE_NORMAL
			case "square" 	=> MODE_SQUARE
			case "tri" 		=> MODE_TRI
			case x => throw new Exception("Unknown mode '%s'".format(x))
		}
	
	
	
	def handleSwathFile(xr:XmlReader, swathFile:File):Unit = {
		var dh = new MzMLDataHandlers(
				setupDataStructures,
				handleSpectrum,
				nc => {},
				c => {})
		
		
		println("HANDLING FILE "+swathFile)
		val binaryFileChannel = 
			if (isIMzML(swathFile)) {
				val binaryFile = toIBD(swathFile)
				println("  BINARY FILE "+binaryFile)
				new FileInputStream(binaryFile).getChannel()
			} else null
		var mzML = MzML.fromFile(xr, dh, binaryFileChannel)
		println("DONE WITH "+swathFile)
		println
		
		val runtime = Runtime.getRuntime
		println(" heap size before writing: "+(runtime.totalMemory / 1000000) + " / "+(runtime.maxMemory / 1000000) + " Mb")
		println()
		
		for (j <- 0 until tramls.length) {
			currTramlOut = j
			var tramlFile = params.tramlFiles(j)
			val outFile = 
				{
				val base = 
					if (tramlFile.toString.toLowerCase.endsWith(".traml"))
						new File(tramlFile.toString.dropRight(6)+".chrom.mzML")
					else
						new File(tramlFile.toString+".chrom.mzML")
				new File(params.outDir, base.getName)
				}
			
			
			println(" writing output file: "+outFile)
			
			mzML.fileDescription.fileContent = CHROM_FILE_CONTENT
			mzML.softwares += SOFTWARE
			mzML.dataProcessings += DATA_PROCESSING
				
			out = XmlWriter(outFile, false)
			var dw = new MzMLDataWriters(
							0,
							ws => Nil,
							transitionSets(j).parts.values.filter(_.times.nocopyResult.nonEmpty).map(_.chroms.length).sum,
							writeChroms
						)
			
			mzML.write(out, dw)
		}
	}
	
	
	
	def setupDataStructures(numSpec:Int) = {
		this.numSpec = numSpec
		
		println("      num spectra: "+numSpec)
		
		if (params.makeBar) {
			println("|                    |")
			print("|")
		}
		
		
    	/*for (i <- 0 until SWATHS_IN_FILE) {
    		times(i) 		= new Array[Double](specPerSwath)
    		specCounters(i) = -1
    	}*/
		ms1Times =
			if (params.subSampleMs1.value == 1) 	
				ChromBuilder()
			else
				AggrChromBuilder(params.subSampleMs1, true)
		
		val l = tramls.length
    	for (j <- 0 until l) {
    		val traml = tramls(j)
    		transitionSets 	+= new TransitionPartitioning(traml.transitions.map(t => new ChromExtract(t)))
    		isotopeSets 	+= traml.transitionGroups.toSeq
	    								.map(t => getIsotopes(traml)(t._2))
	    								.flatten.sortBy(_.q1)
	    								.map(iso => new ChromExtract(iso, params.subSampleMs1))
	    	//isotopeChromSets(j) = isotopeSets(j).map(_ => new Array[Double](specPerSwath))
	    }
	}
	
	
	
	def getIsotopes(traml:GhostTraML)(gts:Seq[GhostTransition]):Seq[Isotope] = {
		if (params.nIsotopes.value == 0)
			return Nil
		
		gts.head.peptide match {
			case None =>
				throw new Exception("Cannot compute isotopes for non-peptide assay! First transition id=" +gts.head.id)
			case Some(pep) =>
				val seq 	= pep.sequence
				try {
					val p 		= UniMod.parseUniModSequence(seq)
					val q1		= gts.head.q1
					val q1z 	= math.round(p.monoisotopicMass() / q1).toDouble
					val id 		= p.getIsotopeDistribution()
					val isotopes = new ArrayBuffer[Isotope]
					for (i <- id.intensities.sorted.takeRight(params.nIsotopes)) {
						val ii = id.intensities.indexOf(i)
						isotopes += new Isotope(seq, q1 + ii / q1z, i)
					}
					isotopes
				} catch {
					case iae:IllegalArgumentException =>
						println("Unable to parse peptide from amino acid sequence '"+seq+"'... ")
						println(iae.getMessage())
						return Nil
				}
		}
			
		
	}
	
	
	
	def handleSpectrum(s:Spectrum):Unit = {
		var gs = GhostSpectrum.fromSpectrum(s)
		
		specCount += 1
		if (params.makeBar) {
			if (barCount < 20 && 20*specCount / numSpec > barCount) {
				print("=")
				barCount += 1
				if (barCount == 20)
					println("|")
			}
		}
		
		if (specCount > numSpec*capFixPercent && !capacitiesFixed) {
			def hint[T](cb:ChromBuilder) = {
				val hl = (cb.nocopyResult.length * numSpec) / (specCount - 2)
				cb.sizeHint(hl)
			}
			
			for (ts <- transitionSets)
				for (partExtract <- ts.parts.values) {
					hint(partExtract.times)
					for (chromExtract <- partExtract.chroms)
						hint(chromExtract.chrom)
				}
			for (set <- isotopeSets)
				for (is <- set)
					hint(is.chrom)
			hint(ms1Times)
			capacitiesFixed = true
		}
		
		val mzs 			= gs.mzs
		val intensities 	= gs.intensities
		val wl = mzs.length
				
		gs.msLevel match {
			case GhostSpectrum.MS1 =>
				ms1SpecCounter += 1
				ms1Times += gs.scanStartTime
				
				for (i <- 0 until isotopeSets.length) {
					val isotopes 	= isotopeSets(i)
					
					var j 	= 0
					var wstart = 0
					var wend = 0
					val jl 	= isotopes.length
			
					while (j < jl) {
						val chromExtract = isotopes(j)
						val q1 = chromExtract.id.q1
						if (j > 0 && q1 < isotopes(j-1).id.q1)
							throw new Exception("isotopes are not sorted in ascending q1 order!")
						var intensity = 0.0
						val cutoff = params.mzThreshold.cutoff(q1)
						
						while (wstart < wl && mzs(wstart) < q1 - cutoff) wstart += 1
						wend = math.max(wstart, wend)
						while (wend < wl && mzs(wend) < q1 + cutoff) wend += 1
						var w = wstart
	
						if (mode == MODE_BEST) {
							var minDiff = Double.MaxValue
							var bestW = -1
							while (w < wend) {
								var d = math.abs(mzs(w) - q1)
								if (d < minDiff) {
									minDiff = d
									bestW = w
								}
								w += 1
							}
							if (minDiff < cutoff)
								intensity = intensities(bestW)
						} else if (mode == MODE_UNIFORM){
							while (w < wend) {
								intensity += intensities(w)
								w += 1
							}
						} else if (mode == MODE_SQUARE) {
							while (w < wend) {
								val d = math.abs(mzs(w) - q1)
								val x = d / cutoff
								intensity += intensities(w) * (1 - (x * x))
								w += 1
							}
						} else if (mode == MODE_NORMAL) {
							while (w < wend) {
								val d = math.abs(mzs(w) - q1)
								val x = d / (2*cutoff)
								intensity += intensities(w) * 2.50662827 * nd.density(x)
								w += 1
							}
						} else if (mode == MODE_TRI) {
							while (w < wend) {
								val d = math.abs(mzs(w) - q1)
								val x = d / cutoff
								intensity += intensities(w) * (1-x)
								w += 1
							}
						}
						
						chromExtract.chrom += intensity
						j += 1
					}
				}
			
			case GhostSpectrum.MS2(precMz, isoWindowOpt) =>
				
				val q1window = getQ1Window(gs)
				
				val t = gs.scanStartTime
		    	for (j <- 0 until tramls.length) {
					val extracts 	= transitionSets(j).getTransitions(q1window)
					extracts.times += t
					
					var i 	= 0
					var wstart = 0
					var wend = 0
					val il 	= extracts.chroms.length
			
					while (i < il) {
						val chromExtract = extracts.chroms(i)
						val q3 = chromExtract.id.q3
						if (i > 0 && q3 < extracts.chroms(i-1).id.q3)
							throw new Exception("transitions not sorted by ascending q3")
						var intensity = 0.0
						val cutoff = params.mzThreshold.cutoff(q3)
						
						while (wstart < wl && mzs(wstart) < q3 - cutoff) wstart += 1
						wend = math.max(wstart, wend)
						while (wend < wl && mzs(wend) < q3 + cutoff) wend += 1
						var w = wstart
	
						if (mode == MODE_BEST) {
							var minDiff = Double.MaxValue
							var bestW = -1
							while (w < wend) {
								val d = math.abs(mzs(w) - q3)
								if (d < minDiff) {
									minDiff = d
									bestW = w
								}
								w += 1
							}
							if (minDiff < cutoff)
								intensity = intensities(bestW)
						} else if (mode == MODE_UNIFORM){
							while (w < wend) {
								intensity += intensities(w)
								w += 1
							}
						} else if (mode == MODE_SQUARE) {
							while (w < wend) {
								val x = math.abs(mzs(w) - q3) / cutoff
								intensity += intensities(w) * (1 - (x * x))
								w += 1
							}
						} else if (mode == MODE_NORMAL) {
							while (w < wend) {
								val d = math.abs(mzs(w) - q3)
								val x = d / (2*cutoff)
								intensity += intensities(w) * 2.50662827 * nd.density(x)
								w += 1
							}
						} else if (mode == MODE_TRI) {
							while (w < wend) {
								val d = math.abs(mzs(w) - q3)
								val x = d / cutoff
								intensity += intensities(w) * (1-x)
								w += 1
							}
						}
						
						chromExtract.chrom += intensity
						i += 1
					}
		    	}
		}
	}
	
	
	def writeFragmentChrom(
			w:XmlWriter,
			index:Int,
			pe:TransitionPartitioning.PartitionExtract, 
			ce:ChromExtract[GhostTransition]
	) = {
		var t = ce.id
		var gc = new GhostChromatogram
		gc.precursor 		= t.q1
		gc.product 			= t.q3
		for (ce <- t.ce)
			gc.collisionEnergy 	= ce
		gc.times 			= pe.times.nocopyResult
		gc.intensities 		= ce.chrom.nocopyResult
		gc.timeDef 			= GhostBinaryDataArray.DataDef(true, false, 
								MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false)
		gc.intensityDef 	= GhostBinaryDataArray.DataDef(true, false, 
								MSNumpress.ACC_NUMPRESS_SLOF, GhostBinaryDataArray.Intensity(), false)
		val chrom 	= gc.toChromatogram(index)
		chrom.id 	= t.id
		chrom.write(w, null)
	}
	
	
	def precursorIsotopeID(ce:ChromExtract[Isotope]) =
		"%s naturally %.2f%s @ %.3f".format(ce.id.seq, ce.id.occurence*100, "%", ce.id.q1)
	
	
	def writeIsotopeChrom(
			w:XmlWriter,
			index:Int,
			times:Seq[Double],
			ce:ChromExtract[Isotope]
	) = {
		def occurenceParam(occ:Double) = {
			val u = new UserParam
			u.name = "isotope occurence"
			u.dataType = Some("xsd:string")
			u.value = Some(occ.toString)
			u
		}
		
		var iso = ce.id
		var gc = new GhostChromatogram
		gc.precursor 		= iso.q1
		gc.product 			= 0.0
		gc.collisionEnergy 	= 0.0
		gc.chromatogram.userParams += occurenceParam(iso.occurence)
		gc.times 			= times
		gc.intensities 		= ce.chrom.nocopyResult
		gc.timeDef 			= GhostBinaryDataArray.DataDef(true, false, 
									MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false)
		gc.intensityDef 	= GhostBinaryDataArray.DataDef(true, false, 
									MSNumpress.ACC_NUMPRESS_SLOF, GhostBinaryDataArray.Intensity(), false)
		val chrom 	= gc.toChromatogram(index)
		chrom.id	= precursorIsotopeID(ce)
		chrom.write(w, null)
	}
	
	
	def writeChroms(w:XmlWriter):Seq[MzML.OffsetRef] = {
		var index = 0
		val tramlChromSet = transitionSets(currTramlOut)
		val fragChromOffsets = 
			for {
				partitionExtract <- tramlChromSet.parts.values
				chromExtract <- partitionExtract.chroms
			} yield {
				val offsetRef = MzML.OffsetRef(chromExtract.id.id, w.byteOffset, None, None)
				writeFragmentChrom(w, index, partitionExtract, chromExtract)
				index += 1
				offsetRef
			}
		
		val precursorIsotopeOffsets =
			if (ms1SpecCounter > 0) {
				val _times = ms1Times.nocopyResult
				val chroms = isotopeSets(currTramlOut)
				for (chromExtract <- chroms) yield {
					val offsetRef = MzML.OffsetRef(precursorIsotopeID(chromExtract), w.byteOffset, None, None)
					writeIsotopeChrom(w, index, _times, chromExtract)
					index += 1
					offsetRef
				}
			} else Nil
		
		fragChromOffsets.toSeq ++ precursorIsotopeOffsets
	}
	
	
	
	
	def gilletQ1Window(gs:GhostSpectrum) = {
		gs.msLevel match {
			case GhostSpectrum.MS1 => throw new Exception("Can't guess swath window for ms1 spectrum!")
			case GhostSpectrum.MS2(precMz, isoWindowOpt) =>
				val swathIndex = math.floor((precMz - 400) / 25).toInt
				TransitionPartitioning.Q1Window(400 + swathIndex*25, 425 + swathIndex*25)
		}
	}
	
	def snoopQ1Window(gs:GhostSpectrum) = {
		gs.msLevel match {
			case GhostSpectrum.MS1 => throw new Exception("Spectrum %d: Can't guess swath window for ms1 spectrum!".format(gs.spectrum.index))
			case GhostSpectrum.MS2(precMz, None) => 
				throw new Exception("Can't snoop window because no isolation window infomation is found!")
			case GhostSpectrum.MS2(precMz, Some(isoWindow)) =>
				TransitionPartitioning.Q1Window(isoWindow.low, isoWindow.high)
		}
	}
	/*
	
	def writeBinaryArray(w:XmlWriter, a:Array[Double], cvParam:CvParam) = {
		var encoded = Base64.codeDoubleArray(a)
		out.startElement(BINARY_DATA_ARRAY)
		out.writeAttribute(ENCODED_LENGTH, encoded.length)
		
		out.startElement(CV_PARAM)
		out.writeAttribute(CV_REF, "MS")
		out.writeAttribute(ACCESSION, BIT_64_ACC)
		out.writeAttribute(NAME, "64-bit float")
		out.endElement
		
		out.startElement(CV_PARAM)
		out.writeAttribute(CV_REF, "MS")
		out.writeAttribute(ACCESSION, NO_COMPRESSION_ACC)
		out.writeAttribute(NAME, "no compression")
		out.endElement
		
		cvParam.write(w)
		
		out.startElement(BINARY)
		out.text(encoded)
		out.endElement
		
		out.endElement
	}
	*/
}
