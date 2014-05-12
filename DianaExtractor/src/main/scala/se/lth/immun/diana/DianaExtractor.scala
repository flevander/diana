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

object DianaExtractor extends CLIApplication {

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
	
	
	
	
	
	var swathFile:File 			= null
	var tramlFiles:Seq[File]	= Nil
	var tramls:Seq[GhostTraML]	= Nil
	var xr:XmlReader 			= null
	var numSpec					= 0
	var specPerSwath			= 0
	var diaWindowMode			= "gillet"
	var getQ1Window:(GhostSpectrum => TransitionPartitioning.Q1Window) = gilletQ1Window _
	var swathWidth 				= 25.0
	var minDiffCutoff			= 0.0
	var minDiffPPM				= 10
	var mode					= "uniform"
	var force					= true
	var noBar					= false
	var nIsotopes				= 0
	var nd 					= new NormalDistribution

	var transitionSets		= new ArrayBuffer[TransitionPartitioning]
	var isotopeSets			= new ArrayBuffer[Seq[ChromExtract[Isotope]]]
	//var isotopeChromSets	= new ArrayBuffer[Array[Array[Double]]]
	var times 				= new Array[Array[Double]](SWATHS_IN_FILE)
	var ms1Times:ChromBuilder = _
	//var specCounters		= new Array[Int](SWATHS_IN_FILE)
	
	var subSampleMs1 		= 1
	var ms1SpecCounter		= -1
	var capacitiesFixed		= false
	var capFixPercent		= 0.1
	
	
	var outFile:File 	= null
	var out:XmlWriter 	= null
	var currTramlOut	= -1
	var outDir:File 	= null
	
	var barCount = 0
	var specCount = 0
	
	def isIMzML(f:File) 	= f.getName.toLowerCase.endsWith(".imzml")
	def toIBD(f:File)		= new File(f.getAbsolutePath().dropRight(5) + "ibd")
	
	
	
	def main(args:Array[String]):Unit = {
		
		var properties = new Properties
    	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
    	
    	
		
		arg("SWATH_MZML", s => {
			swathFile = new File(s)
		})
		
		rest("TRAML_FILES", strs => {
			tramlFiles = strs.map(str => new File(str))
			tramls = tramlFiles.map(f => GhostTraML.fromFile(new XmlReader(
											new BufferedReader(new FileReader(f))
										)))
			}, false)
		
		opt("ms-resolution", 
				"Interpreted as ppm if integer, and as m/z if float", 
				s => {
					try {
						minDiffPPM = s.toInt
					} catch {
						case _:Throwable => {
							minDiffCutoff = s.toDouble
							minDiffPPM = 0
						}
					}
				},
				"X")
		
		opt("mode", 
				"Mode for XIC extraction, best|uniform|normal|square|tri (default: uniform) ", 
				s => 
					if (Array("best", "uniform", "normal", "square", "tri").contains(s))
						mode = s
					else 
						throw new Exception("Unknown extraction mode '"+s+"'"), 
				"X")
		
		opt("dia-windows", 
				"Mode of dia window detection, gillet | snoop(guess from file) (default "+diaWindowMode+" (1 MS1 + 32 MS2))", 
				s => 
					if (Array("gillet", "snoop").contains(s)) {
						diaWindowMode = s
						if (s == "gillet") 		getQ1Window = gilletQ1Window _
						else if (s == "snoop") 	getQ1Window = snoopQ1Window _
					} else 
						throw new Exception("Unknown dia window mode '"+s+"'"), "X")
		
		opt("isotopes", 
				"The number highest precursor isotopes to extract (default 0)", 
				s => nIsotopes = s.toInt, "X")
		
		opt("force", 
				"Ignore missing attributes in mzML files", 
				s => force = true)
		
		opt("no-bar", 
				"don't print process bar", 
				s => noBar = true)
		
		opt("sub-sample-ms1", 
				"Sum k consequent ms1 spectra and store as one value", 
				s => subSampleMs1 = s.toInt, "k")
		
		opt("out-dir", 
				"Directory where chromatogram mzML files should be stored.", 
				s => outDir = new File(s), "X")
		
		val before 		= System.currentTimeMillis
    	val name 		= properties.getProperty("pom.name")
    	val version 	= properties.getProperty("pom.version")
    	
		try {
			parseArgs(name + " "+version, args)
		} catch {
			case cae:CommandlineArgumentException => return
		}
		
		println(name + " "+version)
    	println("  swath mzML file: " + swathFile)
    	println("     num isotopes: " + nIsotopes)
    	println("   sub sample ms1: " + subSampleMs1)
    	println("            force: " + force)
    	println("     max diff PPM: " + minDiffPPM)
    	println("      max diff Da: " + minDiffCutoff)
    	println("  extraction mode: " + mode)
		println("      traML files:")
    	tramlFiles.foreach(f => println("    "+f))
    	println()
    	
		xr = new XmlReader(
				if (swathFile.getName.toLowerCase.endsWith(".gz"))
					new BufferedReader(new InputStreamReader(
						new GZIPInputStream(new FileInputStream(swathFile))))
				else
					new BufferedReader(new FileReader(swathFile))
			)
		xr.force = force
		handleSwathFile(xr, swathFile)
		 
		val after = System.currentTimeMillis
		val nSpectraPerPartition = transitionSets.head.parts.values.map(_.times.nocopyResult.length)
		val runt = Runtime.getRuntime
		println("        heap size: "+(runt.totalMemory / 1000000) + " / "+(runt.maxMemory / 1000000) + " Mb")
		println("    n ms1 spectra: "+ms1SpecCounter)
		println("    n dia windows: "+transitionSets.head.parts.size)
		println("spectra per swath: "+nSpectraPerPartition.min + " - "+nSpectraPerPartition.max)
		println("       time taken: "+niceTiming(after-before))
		println("     spectra read: "+transitionSets.head.parts.values.map(_.times.nocopyResult.length).sum)
	}
	
	
	
	def niceTiming(t:Long) = {
		val ms = t % 1000
		var x = t / 1000
		val s = x % 60
		x = x / 60
		val m = x & 60
		x = x / 60
		val h = x % 24
		val d = x / 24
		"%d days %02d:%02d:%02d.%ds".format(d, h, m, s, ms)
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
		
		val runt = Runtime.getRuntime
		println(" heap size before writing: "+(runt.totalMemory / 1000000) + " / "+(runt.maxMemory / 1000000) + " Mb")
		println()
		
		for (j <- 0 until tramls.length) {
			currTramlOut = j
			var tramlFile = tramlFiles(j)
			val outFile = 
				{
				val base = 
					if (tramlFile.toString.toLowerCase.endsWith(".traml"))
						new File(tramlFile.toString.dropRight(6)+".chrom.mzML")
					else
						new File(tramlFile.toString+".chrom.mzML")
				if (outDir != null)
					new File(outDir, base.getName())
				else
					base
				}
			
			
			println(" writing output file: "+outFile)
			
			mzML.fileDescription.fileContent = CHROM_FILE_CONTENT
			mzML.softwares += SOFTWARE
			mzML.dataProcessings += DATA_PROCESSING
				
			out = new XmlWriter(new BufferedWriter(new FileWriter(outFile)))
			var dw = new MzMLDataWriters(
							0,
							ws => {},
							transitionSets(j).parts.values.filter(_.times.nocopyResult.nonEmpty).map(_.chroms.length).sum,
							writeChroms
						)
			
			mzML.write(out, dw)
		}
	}
	
	
	
	def setupDataStructures(numSpec:Int) = {
		this.numSpec = numSpec
		
		println("      num spectra: "+numSpec)
		
		if (!noBar) {
			println("|                    |")
			print("|")
		}
		
		
    	/*for (i <- 0 until SWATHS_IN_FILE) {
    		times(i) 		= new Array[Double](specPerSwath)
    		specCounters(i) = -1
    	}*/
		ms1Times =
			if (subSampleMs1 == 1) 	ChromBuilder()
			else					AggrChromBuilder(subSampleMs1, true)
		
		val l = tramls.length
    	for (j <- 0 until l) {
    		val traml = tramls(j)
    		transitionSets 	+= new TransitionPartitioning(traml.transitions.map(t => new ChromExtract(t)))
    		isotopeSets 	+= traml.transitionGroups.toSeq
	    								.map(t => getIsotopes(traml)(t._2))
	    								.flatten.sortBy(_.q1)
	    								.map(iso => new ChromExtract(iso, subSampleMs1))
	    	//isotopeChromSets(j) = isotopeSets(j).map(_ => new Array[Double](specPerSwath))
	    }
	}
	
	
	
	def getIsotopes(traml:GhostTraML)(gts:Seq[GhostTransition]):Seq[Isotope] = {
		if (nIsotopes == 0)
			return Nil
		
		val seq 	= traml.peptides(gts.head.peptideRef).sequence
		try {
			val p 		= UniMod.parseUniModSequence(seq)
			val q1		= gts.head.q1
			val q1z 	= math.round(p.monoisotopicMass() / q1).toDouble
			val id 		= p.getIsotopeDistribution()
			val isotopes = new ArrayBuffer[Isotope]
			for (i <- id.intensities.sorted.takeRight(nIsotopes)) {
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
	
	
	
	def handleSpectrum(s:Spectrum):Unit = {
		var gs = GhostSpectrum.fromSpectrum(s)
		
		specCount += 1
		if (!noBar) {
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
				
		if (gs.msLevel == 1) {
			
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
					val ppmCutoff = q1 / 1000000 * minDiffPPM
					var intensity = 0.0
					val cutoff = math.max(ppmCutoff, minDiffCutoff)
					
					while (wstart < wl && mzs(wstart) < q1 - cutoff) wstart += 1
					wend = math.max(wstart, wend)
					while (wend < wl && mzs(wend) < q1 + cutoff) wend += 1
					var w = wstart

					if (mode == "best") {
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
					} else if (mode == "uniform"){
						while (w < wend) {
							intensity += intensities(w)
							w += 1
						}
					} else if (mode == "square") {
						while (w < wend) {
							val d = math.abs(mzs(w) - q1)
							val x = d / cutoff
							intensity += intensities(w) * (1 - (x * x))
							w += 1
						}
					} else if (mode == "normal") {
						while (w < wend) {
							val d = math.abs(mzs(w) - q1)
							val x = d / (2*cutoff)
							intensity += intensities(w) * 2.50662827 * nd.density(x)
							w += 1
						}
					} else if (mode == "tri") {
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
			
		} else if (gs.msLevel == 2) {
			
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
					val ppmCutoff = q3 / 1000000 * minDiffPPM
					var intensity = 0.0
					val cutoff = math.max(ppmCutoff, minDiffCutoff)
					
					while (wstart < wl && mzs(wstart) < q3 - cutoff) wstart += 1
					wend = math.max(wstart, wend)
					while (wend < wl && mzs(wend) < q3 + cutoff) wend += 1
					var w = wstart

					if (mode == "best") {
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
					} else if (mode == "uniform"){
						while (w < wend) {
							intensity += intensities(w)
							w += 1
						}
					} else if (mode == "square") {
						while (w < wend) {
							val x = math.abs(mzs(w) - q3) / cutoff
							intensity += intensities(w) * (1 - (x * x))
							w += 1
						}
					} else if (mode == "normal") {
						while (w < wend) {
							val d = math.abs(mzs(w) - q3)
							val x = d / (2*cutoff)
							intensity += intensities(w) * 2.50662827 * nd.density(x)
							w += 1
						}
					} else if (mode == "tri") {
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
		gc.collisionEnergy 	= t.ce
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
		chrom.id	= "%s naturally %.2f%s @ %.3f".format(ce.id.seq, ce.id.occurence*100, "%", ce.id.q1)
		chrom.write(w, null)
	}
	
	
	def writeChroms(w:XmlWriter) = {
		var index = 0
		val tramlChromSet = transitionSets(currTramlOut)
		for (partitionExtract <- tramlChromSet.parts.values)
			for (chromExtract <- partitionExtract.chroms) {
				writeFragmentChrom(w, index, partitionExtract, chromExtract)
				index += 1
			}
		
		val _times = ms1Times.nocopyResult
		val chroms = isotopeSets(currTramlOut)
		if (ms1SpecCounter > 0) {
			for (chromExtract <- chroms) {
				writeIsotopeChrom(w, index, _times, chromExtract)
				index += 1
			}
		}
	}
	
	
	
	
	def gilletQ1Window(gs:GhostSpectrum) = {
		val swathIndex = math.floor((gs.q1 - 400) / 25).toInt
		TransitionPartitioning.Q1Window(400 + swathIndex*25, 425 + swathIndex*25)
	}
	
	def snoopQ1Window(gs:GhostSpectrum) = {
		val cvs = gs.spectrum.precursors.head.isolationWindow.get.cvParams
		val isolationWindowTargetMz = cvs.find(_.accession == "MS:1000827").get.value.get.toDouble
		val isolationWindowLowerOffset = cvs.find(_.accession == "MS:1000828").get.value.get.toDouble
		val isolationWindowUpperOffset = cvs.find(_.accession == "MS:1000829").get.value.get.toDouble
		TransitionPartitioning.Q1Window(
				isolationWindowTargetMz - isolationWindowLowerOffset, 
				isolationWindowTargetMz + isolationWindowUpperOffset)
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
