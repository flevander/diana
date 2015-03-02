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
import se.lth.immun.collection.numpress.NumLinArray
import se.lth.immun.collection.numpress.NumSlofArray

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.ArrayBuilder

import org.apache.commons.math3.distribution.NormalDistribution

object Tracer extends CLIApplication {

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
			val id:String,
			val q1:Double,
			val occurence:Double
	) {}
	
	case class FragIndex(traml:Int, channel:Int, id:GhostTransition) {
		override def toString = "%s|%.2f-%.2f".format(channel, id.q1, id.q3)
	}
	
	case class ChromGroupIndex(traml:Int, channelGroup:Int, q1:Double) {
		override def toString = "%s|%.2f".format(channelGroup, q1)
	}
	
	
	
	
	
	var mzMLFile:File 			= null
	var tramlFiles:Seq[File]	= Nil
	var tramls:Seq[GhostTraML]	= Nil
	var xr:XmlReader 			= null
	var numSpec					= 0
	var specPerSwath			= 0
	var diaWindowMode			= "snoop"
	var getQ1Windows:(GhostSpectrum => Seq[ChannelLookup.Range]) = snoopQ1Windows _
	var swathWidth 				= 25.0
	var minDiffCutoff			= 0.0
	var minDiffPPM				= 10
	var mode					= "uniform"
	var force					= true
	var noBar					= false
	var nd 					= new NormalDistribution

	var transitionSets		= new ArrayBuffer[Seq[ChannelGroupTrace[Double, GhostTransition]]]
	var isotopeSets			= new ArrayBuffer[Seq[ChannelTrace[Isotope]]]
	var ms1Times:AggrBuilder = _
	//var specCounters		= new Array[Int](SWATHS_IN_FILE)
	
	var chromGroupLookup:ChannelLookup[ChromGroupIndex] = _
	var allIsos:Seq[ChannelTrace[Isotope]] = _
	
	var subSampleMs1 		= 1
	var capacitiesFixed		= false
	var capFixPercent		= 0.1
	
	
	var outFile:File 	= null
	var out:XmlWriter 	= null
	var currTramlOut	= -1
	var outDir:File 	= null
	
	var ms1SpecCounter		= 0
	var ms2SpecCounter		= 0
	var specCount = 0
	var barCount = 0
	
	def isIMzML(f:File) 	= f.getName.toLowerCase.endsWith(".imzml")
	def toIBD(f:File)		= new File(f.getAbsolutePath().dropRight(5) + "ibd")
	
	
	
	def main(args:Array[String]):Unit = {
		
		var properties = new Properties
    	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
    	
    	
		
		arg("MZML", s => {
			mzMLFile = new File(s)
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
				"Mode of dia window detection, gillet | snoop(guess from file) (default "+diaWindowMode+")", 
				s => 
					if (Array("gillet", "snoop").contains(s)) {
						diaWindowMode = s
						if (s == "gillet") 		getQ1Windows = gilletQ1Window _
						else if (s == "snoop") 	getQ1Windows = snoopQ1Windows _
					} else 
						throw new Exception("Unknown dia window mode '"+s+"'"), "X")
		
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
    	println("  swath mzML file: " + mzMLFile)
    	println("   sub sample ms1: " + subSampleMs1)
    	println("            force: " + force)
    	println("     max diff PPM: " + minDiffPPM)
    	println("      max diff Da: " + minDiffCutoff)
    	println("  extraction mode: " + mode)
		println("      traML files:")
    	tramlFiles.foreach(f => println("    "+f))
    	println()
    	
		xr = new XmlReader(
				if (mzMLFile.getName.toLowerCase.endsWith(".gz"))
					new BufferedReader(new InputStreamReader(
						new GZIPInputStream(new FileInputStream(mzMLFile))))
				else
					new BufferedReader(new FileReader(mzMLFile))
			)
		xr.force = force
		handleSwathFile(xr, mzMLFile)
		 
		val after = System.currentTimeMillis
		val runt = Runtime.getRuntime
		println("        heap size: "+(runt.totalMemory / 1000000) + " / "+(runt.maxMemory / 1000000) + " Mb")
		println("    n ms1 spectra: "+ms1SpecCounter)
		println("       time taken: "+niceTiming(after-before))
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
		"%d days %02d:%02d:%02d.%03ds".format(d, h, m, s, ms)
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
							transitionSets(j).map(_.subChannels.size).sum + isotopeSets(j).length,
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
			new AggrBuilder(new NumLinArray(10000.0, 20), subSampleMs1, true)
		
		val l = tramls.length
    	for (j <- 0 until l) {
    		val traml = tramls(j)
    		transitionSets 	+= 
    			(for {
	    			(q1, transes) <- traml.transitions.groupBy(_.q1).toSeq
	    		} yield new ChannelGroupTrace(q1, transes))
    		isotopeSets 	+= traml.includes.map(gt => new ChannelTrace(toIso(gt), subSampleMs1))
	    	//isotopeChromSets(j) = isotopeSets(j).map(_ => new Array[Double](specPerSwath))
	    }
		
		val chromGroupIndices = 
			for {
				(tset, i) <- transitionSets.zipWithIndex
				(cgt, j) <- tset.zipWithIndex
			} yield ChromGroupIndex(i, j, cgt.id)
		
		chromGroupLookup = new ChannelLookup(chromGroupIndices, _.q1)
		allIsos = isotopeSets.flatten.sortBy(_.id.q1)
	}
	
	
	
	def toIso(gt:GhostTarget):Isotope = {
		new Isotope(gt.id, gt.q1, gt.intensity)
	}
	
	
	/*
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
	*/
	
	
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
		/*
		if (specCount > numSpec*capFixPercent && !capacitiesFixed) {
			def hint[T](cb:TraceBuilder) = {
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
		*/
		val mzs 			= gs.mzs
		val intensities 	= gs.intensities
		val wl = mzs.length
				
		if (gs.msLevel == 1) {
			
			ms1SpecCounter += 1
			ms1Times += gs.scanStartTime
			
			var j 	= 0
			var wstart = 0
			var wend = 0
			val jl 	= allIsos.length
		
			while (j < jl) {
				val channelTrace = allIsos(j)
				val q1 = channelTrace.id.q1
				if (j > 0 && q1 < allIsos(j-1).id.q1)
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
				
				channelTrace.intensity += intensity
				j += 1
			}
			
		} else if (gs.msLevel == 2) {
			
			ms2SpecCounter += 1
			
			val q1windows 	= getQ1Windows(gs)
			val t 			= gs.scanStartTime
	    	val cgIndices	= chromGroupLookup.fromRanges(q1windows).toSeq
	    	val traces		= cgIndices.flatMap(cgi => 
	    						transitionSets(cgi.traml)(cgi.channelGroup).subChannels
	    					).sortBy(_._1.q3)
			
	    	//println(q1windows.mkString(" "))
	    	//println(fragIndices.mkString(" "))
	    	
			var i 		= 0
			var wstart 	= 0
			var wend 	= 0
			val il 		= traces.length
	
			for (cgi <- cgIndices) {
				transitionSets(cgi.traml)(cgi.channelGroup).time += t
			}
			
			while (i < il) {
				val (gt, intensityTrace) = traces(i)
				val q3 = gt.q3
				if (i > 0 && q3 < traces(i-1)._1.q3)
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
							bestW 	= w
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
				
				intensityTrace += intensity
				i += 1
			}
		}
	}
	
	
	def writeFragmentChrom(
			w:XmlWriter,
			index:Int,
			gt:GhostTransition,
			times:NumLinArray,
			intensities:NumSlofArray
	) = {
		val gc = new GhostChromatogram
		gc.precursor 		= gt.q1
		gc.product 			= gt.q3
		gc.collisionEnergy 	= gt.ce
		
		val c = new Chromatogram
		c.id = gt.id
		c.index = index
		c.defaultArrayLength = times.length
		c.precursor = gc.chromatogram.precursor
		c.product 	= gc.chromatogram.product
		
		c.binaryDataArrays += ToMzML.toBinaryDataArray(
									times.ba, 
									GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false))
		c.binaryDataArrays += ToMzML.toBinaryDataArray(
									intensities.ba,
									GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_SLOF, GhostBinaryDataArray.Intensity(), false))
		c.write(w, null)
	}
	
	
	def writeIsotopeChrom(
			w:XmlWriter,
			index:Int,
			times:AggrBuilder,
			ct:ChannelTrace[Isotope]
	) = {
		def occurenceParam(occ:Double) = {
			val u = new UserParam
			u.name = "isotope occurence"
			u.dataType = Some("xsd:string")
			u.value = Some(occ.toString)
			u
		}
		
		val iso = ct.id
		val gc = new GhostChromatogram
		gc.precursor 		= iso.q1
		gc.product 			= 0.0
		gc.collisionEnergy 	= 0.0
		
		val c = new Chromatogram
		c.id = iso.id
		c.index = index
		c.defaultArrayLength = times.a.length
		c.userParams += occurenceParam(iso.occurence)
		c.precursor = gc.chromatogram.precursor
		c.product 	= gc.chromatogram.product
		
		c.binaryDataArrays += ToMzML.toBinaryDataArray(
									times.a.ba, 
									GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false))
		c.binaryDataArrays += ToMzML.toBinaryDataArray(
									ct.intensity.a.ba,
									GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_SLOF, GhostBinaryDataArray.Intensity(), false))
		c.write(w, null)
	}
	
	
	def writeChroms(w:XmlWriter) = {
		var index = 0
		val tramlChannelGroups = transitionSets(currTramlOut)
		for {
			cgt <- tramlChannelGroups
			(gt, trace) <- cgt.subChannels
		} {
			writeFragmentChrom(w, index, gt, cgt.time, trace)
			index += 1
		}
		
		val tramlIsotopeTraces = isotopeSets(currTramlOut)
		if (ms1SpecCounter > 0) {
			for (trace <- tramlIsotopeTraces) {
				writeIsotopeChrom(w, index, ms1Times, trace)
				index += 1
			}
		}
	}
	
	
	
	
	def gilletQ1Window(gs:GhostSpectrum) = {
		val swathIndex = math.floor((gs.q1 - 400) / 25).toInt
		List(ChannelLookup.Range(400 + swathIndex*25, 425 + swathIndex*25))
	}
	
	def snoopQ1Windows(gs:GhostSpectrum) = {
		gs.spectrum.precursors.map(pc => {
			val cvs = pc.isolationWindow.get.cvParams
			val isolationWindowTargetMz = cvs.find(_.accession == "MS:1000827").get.value.get.toDouble
			val isolationWindowLowerOffset = cvs.find(_.accession == "MS:1000828").get.value.get.toDouble
			val isolationWindowUpperOffset = cvs.find(_.accession == "MS:1000829").get.value.get.toDouble
			ChannelLookup.Range(
				isolationWindowTargetMz - isolationWindowLowerOffset, 
				isolationWindowTargetMz + isolationWindowUpperOffset)
		})
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
