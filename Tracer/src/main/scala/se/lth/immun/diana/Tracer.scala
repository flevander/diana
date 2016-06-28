package se.lth.immun.diana

import se.lth.immun.xml.XmlReader
import se.lth.immun.xml.XmlWriter
import se.jt.CLIApp
import se.jt.CLIBar

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

object Tracer extends CLIApp {

	import MzML._
	
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
		id = "Tracer"
		version = "0.3.3"
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
		
		id = "DIAExtraction"
		processingMethods += pm
	}
	
	class Isotope(
			val id:String,
			val q1:Double,
			val occurence:Double
	) {}
	
	//case class FragIndex(traml:Int, channel:Int, id:GhostTransition) {
	//	override def toString = "%s|%.2f-%.2f".format(channel, id.q1, id.q3)
	//}
	
	case class ChromGroupIndex(traml:Int, channelGroup:Int, q1:Double) {
		override def toString = "%s|%.2f".format(channelGroup, q1)
	}
	
	
	
	
	
	var mzMLFile:File 			= null
	var tramlFiles:Seq[File]	= Nil
	var tramls:Seq[GhostTraML]	= Nil
	var numSpec					= 0
	var getQ1Windows:(GhostSpectrum => Seq[ChannelLookup.Range]) = snoopQ1Windows _
	
	
	var transitionSets		= new ArrayBuffer[Seq[ChannelGroupTrace[Double, GhostTransition]]]
	var isotopeSets			= new ArrayBuffer[Seq[ChannelTrace[Isotope]]]
	var ms1Times:AggrBuilder = _
	
	var chromGroupLookup:ChannelLookup[ChromGroupIndex] = _
	var allIsos:Seq[ChannelTrace[Isotope]] = _
	
	var nd 					= new NormalDistribution
	
	var params:TracerParams = _
	
	var ms1SpecCounter		= 0
	var ms2SpecCounter		= 0
	var specCount = 0
	var cliBar = new CLIBar
	
	def isIMzML(f:File) 	= f.getName.toLowerCase.endsWith(".imzml")
	def toIBD(f:File)		= new File(f.getAbsolutePath().dropRight(5) + "ibd")
	
	var analysisTime = 0L
	var writeTime = 0L
	
	
	def main(args:Array[String]):Unit = {
		
		var properties = new Properties
    	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
    	val before 		= System.currentTimeMillis
    	val name 		= properties.getProperty("pom.name")
    	val version 	= properties.getProperty("pom.version")
    	
    	params = new TracerParams(name, version)
    	
		failOnError(parseArgs(name, version, args, params, List("mzML"), Some("tramls")))
    	failOnError(params.parseMinDiff)
    	
    	mzMLFile = new File(params.mzML)
		tramlFiles = params.tramls.map(t => new File(t))
    	tramls = tramlFiles.map(f => GhostTraML.fromFile(new XmlReader(
											new BufferedReader(new FileReader(f))
										)))
		
		println(name + " "+version)
    	println("    dia mzML file: " + mzMLFile)
    	println("   sub sample ms1: " + params.subSampleMs1.value)
    	println("            force: " + params.force.value)
    	println("     max diff PPM: " + params.minDiffPPM)
    	println("      max diff Da: " + params.minDiffCutoff)
    	println("  extraction mode: " + params.mode.value)
		println("      traML files:")
    	tramlFiles.foreach(f => println("    "+f))
    	println()
    	
		val xr = new XmlReader(
				if (mzMLFile.getName.toLowerCase.endsWith(".gz"))
					new BufferedReader(new InputStreamReader(
						new GZIPInputStream(new FileInputStream(mzMLFile))))
				else
					new BufferedReader(new FileReader(mzMLFile))
			)
		xr.force = params.force
		handleSwathFile(xr, mzMLFile)
		 
		val after = System.currentTimeMillis
		val runt = Runtime.getRuntime
		println("        heap size: "+(runt.totalMemory / 1000000) + " / "+(runt.maxMemory / 1000000) + " Mb")
		println("    n ms1 spectra: "+ms1SpecCounter)
		println(" total time taken: "+niceTiming(after-before))
		if (params.profiling) {
			println("    analysis time: "+niceTiming(analysisTime))
			println("      write taken: "+niceTiming(writeTime))
		}
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
			
		val t0 = System.currentTimeMillis
		var mzML = MzML.fromFile(xr, dh, binaryFileChannel)
		val t1 = System.currentTimeMillis
		
		val runt = Runtime.getRuntime
		println
		println(" heap size before writing: "+(runt.totalMemory / 1000000) + " / "+(runt.maxMemory / 1000000) + " Mb")
		println
		
		for (j <- 0 until tramls.length) {
			var tramlFile = tramlFiles(j)
			val outFile = 
				{
				val base = 
					if (tramlFile.toString.toLowerCase.endsWith(".traml"))
						new File(tramlFile.toString.dropRight(6)+".chrom.mzML")
					else
						new File(tramlFile.toString+".chrom.mzML")
				if (params.outDir.value != "-")
					new File(params.outDir, base.getName())
				else
					base
				}
			
			
			println(" writing output file: "+outFile)
			
			mzML.fileDescription.fileContent = CHROM_FILE_CONTENT
			mzML.softwares += SOFTWARE
			mzML.dataProcessings += DATA_PROCESSING
			
			val out = XmlWriter.apply(outFile, false, 4096)
			var dw = new MzMLDataWriters(
							0,
							ws => {Nil},
							transitionSets(j).map(_.subChannels.size).sum + isotopeSets(j).length,
							writeChroms(j)
						)
			
			mzML.write(out, dw)
		}
		val t2 = System.currentTimeMillis
		analysisTime = t1 - t0
		writeTime = t2 - t1
	}
	
	
	
	def setupDataStructures(numSpec:Int) = {
		this.numSpec = numSpec
		
		println("      num spectra: "+numSpec)
		
		if (!params.verbose) {
			println(cliBar.reference)
			print(cliBar.update(0))
		}
		
		ms1Times = new AggrBuilder(new NumLinArray(params.timeFixedPoint, 20), params.subSampleMs1, false)
		//	if (params.subSampleMs1.value == 1) new ArrayBuffer[Double]
		//	else new AggrArrayBuffer(params.subSampleMs1, true)
		
		val l = tramls.length
    	for (j <- 0 until l) {
    		val traml = tramls(j)
    		transitionSets 	+= 
    			(for {
	    			(q1, transes) <- traml.transitions.groupBy(_.q1).toSeq
	    		} yield new ChannelGroupTrace(q1, transes, params.timeFixedPoint, params.intFixedPoint))
    		isotopeSets 	+= traml.includes.map(gt => new ChannelTrace(toIso(gt), params.subSampleMs1, false, params.intFixedPoint))
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
	
	
	def handleSpectrum(s:Spectrum):Unit = {
		var gs = GhostSpectrum.fromSpectrum(s)
		
		import TracerParams._
		
		specCount += 1
		if (!params.verbose) 
			print(cliBar.update(specCount, numSpec))
		
		val mzs 			= gs.mzs
		val intensities 	= gs.intensities
		val wl = mzs.length
		val mode = params.imode
				
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
				val ppmCutoff = q1 / 1000000 * params.minDiffPPM
				var intensity = 0.0
				val cutoff = math.max(ppmCutoff, params.minDiffCutoff)
				
				while (wstart < wl && mzs(wstart) < q1 - cutoff) wstart += 1
				wend = math.max(wstart, wend)
				while (wend < wl && mzs(wend) < q1 + cutoff) wend += 1
				var w = wstart

				
				if (mode == BEST) {
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
				} else if (mode == UNIFORM){
					while (w < wend) {
						intensity += intensities(w)
						w += 1
					}
				} else if (mode == SQUARE) {
					while (w < wend) {
						val d = math.abs(mzs(w) - q1)
						val x = d / cutoff
						intensity += intensities(w) * (1 - (x * x))
						w += 1
					}
				} else if (mode == NORMAL) {
					while (w < wend) {
						val d = math.abs(mzs(w) - q1)
						val x = d / (2*cutoff)
						intensity += intensities(w) * 2.50662827 * nd.density(x)
						w += 1
					}
				} else if (mode == TRI) {
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
			
			val t 			= gs.scanStartTime
			val q1windows 	= getQ1Windows(gs)
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
				val ppmCutoff = q3 / 1000000 * params.minDiffPPM
				var intensity = 0.0
				val cutoff = math.max(ppmCutoff, params.minDiffCutoff)
				
				while (wstart < wl && mzs(wstart) < q3 - cutoff) wstart += 1
				wend = math.max(wstart, wend)
				while (wend < wl && mzs(wend) < q3 + cutoff) wend += 1
				var w = wstart

				if (mode == BEST) {
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
				} else if (mode == UNIFORM){
					while (w < wend) {
						intensity += intensities(w)
						w += 1
					}
				} else if (mode == SQUARE) {
					while (w < wend) {
						val x = math.abs(mzs(w) - q3) / cutoff
						intensity += intensities(w) * (1 - (x * x))
						w += 1
					}
				} else if (mode == NORMAL) {
					while (w < wend) {
						val d = math.abs(mzs(w) - q3)
						val x = d / (2*cutoff)
						intensity += intensities(w) * 2.50662827 * nd.density(x)
						w += 1
					}
				} else if (mode == TRI) {
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
									times.bytes, 
									GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false))
		c.binaryDataArrays += ToMzML.toBinaryDataArray(
									intensities.bytes,
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
									times.a.bytes, 
									GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false))
		c.binaryDataArrays += ToMzML.toBinaryDataArray(
									ct.intensity.a.bytes,
									GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_SLOF, GhostBinaryDataArray.Intensity(), false))
		c.write(w, null)
	}
	/*
	
	*/
	
	
	/*
	def writeFragmentChrom(
			w:XmlWriter,
			index:Int,
			gt:GhostTransition,
			times:Seq[Double],
			intensities:Seq[Double]
	) = {
		var gc = new GhostChromatogram
		gc.precursor 		= gt.q1
		gc.product 			= gt.q3
		gc.collisionEnergy 	= gt.ce
		gc.times 			= times
		gc.intensities 		= intensities
		gc.timeDef 			= GhostBinaryDataArray.DataDef(true, false, 
								MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false)
		gc.intensityDef 	= GhostBinaryDataArray.DataDef(true, false, 
								MSNumpress.ACC_NUMPRESS_SLOF, GhostBinaryDataArray.Intensity(), false)
		val chrom 	= gc.toChromatogram(index)
		chrom.id 	= gt.id
		chrom.write(w, null)
	}
	
	
	def writeIsotopeChrom(
			w:XmlWriter,
			index:Int,
			times:Seq[Double],
			ch:ChannelTrace[Isotope]
	) = {
		def occurenceParam(occ:Double) = {
			val u = new UserParam
			u.name = "isotope occurence"
			u.dataType = Some("xsd:string")
			u.value = Some(occ.toString)
			u
		}
		
		var iso = ch.id
		var gc = new GhostChromatogram
		gc.precursor 		= iso.q1
		gc.product 			= 0.0
		gc.collisionEnergy 	= 0.0
		gc.chromatogram.userParams += occurenceParam(iso.occurence)
		gc.times 			= times
		gc.intensities 		= ch.intensity
		gc.timeDef 			= GhostBinaryDataArray.DataDef(true, false, 
									MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false)
		gc.intensityDef 	= GhostBinaryDataArray.DataDef(true, false, 
									MSNumpress.ACC_NUMPRESS_SLOF, GhostBinaryDataArray.Intensity(), false)
		val chrom 	= gc.toChromatogram(index)
		chrom.id	= iso.id
		chrom.write(w, null)
	}
	*/
	
	
	
	def writeChroms(currTramlOut:Int)(w:XmlWriter):Seq[OffsetRef] = {
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
		return Nil
	}
	
	
	
	/**
	def gilletQ1Window(gs:GhostSpectrum) = {
		val swathIndex = math.floor((gs.q1 - 400) / 25).toInt
		List(ChannelLookup.Range(400 + swathIndex*25, 425 + swathIndex*25))
	}*/
	
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
}
