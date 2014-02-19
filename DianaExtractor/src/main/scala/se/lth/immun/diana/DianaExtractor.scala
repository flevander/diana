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
import java.util.Calendar
import java.util.Properties
import java.text.SimpleDateFormat

import se.lth.immun.mzml._
import se.lth.immun.mzml.ghost._
import se.lth.immun.traml.ghost._
import se.lth.immun.chem._

import ms.numpress.MSNumpress

import scala.collection.mutable.ArrayBuffer

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
			val q1:Double,
			val occurence:Double
	) {}
	
	
	
	
	
	var swathFile:File 			= null
	var tramlFiles:Seq[File]	= Nil
	var tramls:Seq[GhostTraML]	= Nil
	var xr:XmlReader 			= null
	var numSpec					= 0
	var specPerSwath			= 0
	var typesOfSpectra			= SWATHS_IN_FILE+1
	var swathWidth 				= 25.0
	var minDiffCutoff			= 0.0
	var minDiffPPM				= 10
	var mode					= "uniform"
	var force					= true
	var nIsotopes				= 0
	var nd 					= new NormalDistribution

	var transitionSets		= new ArrayBuffer[Array[Array[GhostTransition]]]
	var chromSets			= new ArrayBuffer[Array[Array[Array[Double]]]]
	var isotopeSets			= new ArrayBuffer[Array[Isotope]]
	var isotopeChromSets	= new ArrayBuffer[Array[Array[Double]]]
	var times 				= new Array[Array[Double]](SWATHS_IN_FILE)
	var ms1Times			= new ArrayBuffer[Double]
	var specCounters		= new Array[Int](SWATHS_IN_FILE)
	var ms1SpecCounter		= -1
	
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
						case _ => {
							minDiffCutoff = s.toDouble
							minDiffPPM = 0
						}
					}
				},
				"X")
		
		opt("mode", 
				"Mode for XIC extraction, best|uniform|normal|square|tri (default: uniform) ", 
				s => 
					if (Array("best", "uniform", "normal", "square", "tri").contains(mode))
						mode = s
					else 
						throw new Exception("Unknown extraction mode '"+s+"'"), 
				"X")
		
		opt("types-of-spectra", 
				"The number of different precursor swaths (default 32+1)", 
				s => typesOfSpectra = s.toInt, "X")
		
		opt("isotopes", 
				"The number highest precursor isotopes to extract (default 0)", 
				s => nIsotopes = s.toInt, "X")
		
		opt("force", 
				"Ignore missing attributes in mzML files", 
				s => force = true)
		
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
    	println("            force: " + force)
    	println("     min diff PPM: " + minDiffPPM)
    	println("      min diff Da: " + minDiffCutoff)
    	println("  extraction mode: " + mode)
		println("      traML files:")
    	tramlFiles.foreach(f => println("    "+f))
    	println()
    	
		xr = new XmlReader(new BufferedReader(new FileReader(swathFile)))
		xr.force = force
		handleSwathFile(xr, swathFile)
		
		val after = System.currentTimeMillis
		println("  time taken: "+(after-before)+"ms")
		println("spectra read: "+specCounters.sum)
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
							chromSets.zip(specCounters).filter(_._2 > -1).map(_._1.length).sum,
							writeChroms
						)
			
			mzML.write(out, dw)
		}
	}
	
	
	
	def setupDataStructures(numSpec:Int) = {
		this.numSpec = numSpec
		specPerSwath = numSpec / typesOfSpectra
		
		println("      num spectra: "+numSpec)
		println(" types of spectra: "+typesOfSpectra)
		println("spectra per swath: "+specPerSwath)
		println("|                    |")
		print("|")

    	for (i <- 0 until SWATHS_IN_FILE) {
    		times(i) 		= new Array[Double](specPerSwath)
    		specCounters(i) = -1
    	}
		
		val l = tramls.length
    	for (j <- 0 until l) {
    		val traml = tramls(j)
    		transitionSets 	+= new Array(SWATHS_IN_FILE)
    		chromSets		+= new Array(SWATHS_IN_FILE)
    		isotopeChromSets += new Array(SWATHS_IN_FILE)
	    	isotopeSets 	+= traml.transitionGroups.toSeq
	    								.map(t => getIsotopes(traml)(t._2))
	    								.flatten.sortBy(_.q1).toArray
	    	isotopeChromSets(j) = isotopeSets(j).map(_ => new Array[Double](specPerSwath))
    		for (i <- 0 until SWATHS_IN_FILE) {
    			var mz = 400 + swathWidth*i
	    		transitionSets(j)(i) = traml.transitions.filter(tr => 
							tr.q1 >= mz && tr.q1 < mz + swathWidth
	    				).toArray
	    		chromSets(j)(i) 		= transitionSets(j)(i).map(_ => new Array[Double](specPerSwath))
    		}
	    }
	}
	
	
	
	def getIsotopes(traml:GhostTraML)(gts:Seq[GhostTransition]):Seq[Isotope] = {
		val seq 	= traml.peptides(gts.head.peptideRef).sequence
		val p 		= new Peptide(seq.map(c => StandardAminoAcid.fromChar(c)).toArray)
		val q1		= gts.head.q1
		val q1z 	= math.round(p.monoisotopicMass() / q1).toDouble
		val id = p.getIsotopeDistribution()
		val isotopes = new ArrayBuffer[Isotope]
		for (i <- id.intensities.sorted.takeRight(nIsotopes)) {
			val ii = id.intensities.indexOf(i)
			isotopes += new Isotope(q1 + ii / q1z, i)
		}
		isotopes
	}
	
	
	
	def handleSpectrum(s:Spectrum):Unit = {
		var gs = GhostSpectrum.fromSpectrum(s)
		
		specCount += 1
		if (barCount < 20 && 20*specCount / numSpec > barCount) {
			print("=")
			barCount += 1
			if (barCount == 20)
				println("|")
		}
		
		val mzs 			= gs.mzs
		val intensities 	= gs.intensities
				
		if (gs.msLevel == 1) {
			
			ms1SpecCounter += 1
			ms1Times += gs.scanStartTime
			for (i <- 0 until isotopeSets.length) {
				val isotopes 	= isotopeSets(i)
				val chroms 		= isotopeChromSets(i)
				
				var j 	= 0
				val jl 	= isotopes.length
		
				while (j < jl) {
					var k = 0
					var kl = mzs.length
					var q1 = isotopes(j).q1
					var ppmCutoff = q1 / 1000000 * minDiffPPM
					var intensity = 0.0
					var cutoff = math.max(ppmCutoff, minDiffCutoff)
					
					if (mode == "best") {
						var minDiff = Double.MaxValue
						var bestK = -1
						while (k < kl) {
							var d = math.abs(mzs(k) - q1)
							if (d < minDiff) {
								minDiff = d
								bestK = j
							}
							k += 1
						}
						if (minDiff < cutoff)
							intensity = intensities(bestK)
					} else if (mode == "uniform"){
						while (k < kl) {
							var d = math.abs(mzs(k) - q1)
							if (d < cutoff)
								intensity += intensities(k)
							k += 1
						}
					} else if (mode == "square") {
						while (k < kl) {
							var d = math.abs(mzs(k) - q1)
							if (d < cutoff) {
								val x = d / cutoff
								intensity += intensities(k) * (1 - (x * x))
							}
							k += 1
						}
					} else if (mode == "normal") {
						while (k < kl) {
							var d = math.abs(mzs(k) - q1)
							if (d < cutoff) {
								val x = d / (2*cutoff)
								intensity += intensities(k) * 2.50662827 * nd.density(x)
							}
							k += 1
						}
					} else if (mode == "tri") {
						while (k < kl) {
							var d = math.abs(mzs(k) - q1)
							if (d < cutoff) {
								val x = d / cutoff
								intensity += intensities(k) * (1-x)
							}
							k += 1
						}
					}
					
					chroms(j)(ms1SpecCounter) = intensity
					j += 1
				}
			}
			
		} else if (gs.msLevel == 2) {
			
			val si = math.floor((gs.q1 - 400) / 25).toInt
			specCounters(si) += 1
			val sc = specCounters(si)
			
			times(si)(sc) 		= gs.scanStartTime
	    	for (j <- 0 until tramls.length) {
				val transitions 	= transitionSets(j)(si)
				val chroms 			= chromSets(j)(si)
				
				var i 	= 0
				val il 	= transitions.length
		
				while (i < il) {
					var k = 0
					var kl = mzs.length
					var q3 = transitions(i).q3
					var ppmCutoff = q3 / 1000000 * minDiffPPM
					var intensity = 0.0
					var cutoff = math.max(ppmCutoff, minDiffCutoff)
					
					if (mode == "best") {
						var minDiff = Double.MaxValue
						var bestK = -1
						while (k < kl) {
							var d = math.abs(mzs(k) - q3)
							if (d < minDiff) {
								minDiff = d
								bestK = j
							}
							k += 1
						}
						if (minDiff < cutoff)
							intensity = intensities(bestK)
					} else if (mode == "uniform"){
						while (k < kl) {
							var d = math.abs(mzs(k) - q3)
							if (d < cutoff)
								intensity += intensities(k)
							k += 1
						}
					} else if (mode == "square") {
						while (k < kl) {
							var d = math.abs(mzs(k) - q3)
							if (d < cutoff) {
								val x = d / cutoff
								intensity += intensities(k) * (1 - (x * x))
							}
							k += 1
						}
					} else if (mode == "normal") {
						while (k < kl) {
							var d = math.abs(mzs(k) - q3)
							if (d < cutoff) {
								val x = d / (2*cutoff)
								intensity += intensities(k) * 2.50662827 * nd.density(x)
							}
							k += 1
						}
					} else if (mode == "tri") {
						while (k < kl) {
							var d = math.abs(mzs(k) - q3)
							if (d < cutoff) {
								val x = d / cutoff
								intensity += intensities(k) * (1-x)
							}
							k += 1
						}
					}
					
					chroms(i)(sc) = intensity
					i += 1
				}
	    	}
		}
	}
	
	
	def writeChroms(w:XmlWriter) = {
		var index = 0
		val tramlChromSet = chromSets(currTramlOut)
		for (iset <- 0 until tramlChromSet.length) {
			if (specCounters(iset) > 0) {
				var chroms = tramlChromSet(iset)
				for (ichrom <- 0 until chroms.length) {
					var t = transitionSets(currTramlOut)(iset)(ichrom)
					var gc = new GhostChromatogram
					gc.precursor 		= t.q1
					gc.product 			= t.q3
					gc.collisionEnergy 	= t.ce
					gc.times 			= times(iset)
					gc.intensities 		= chroms(ichrom)
					gc.timeDef 			= GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false)
					gc.intensityDef 	= GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_SLOF, GhostBinaryDataArray.Intensity(), false)
					gc.toChromatogram(index).write(w, null)
					index += 1
				}
			}
		}
		
		def occurenceParam(occ:Double) = {
			val u = new UserParam
			u.name = "isotope occurence"
			u.dataType = Some("xsd:string")
			u.value = Some(occ.toString)
			u
		}
		
		val _times = ms1Times.toArray
		val chroms = isotopeChromSets(currTramlOut)
		if (ms1SpecCounter > 0) {
			for (ichrom <- 0 until chroms.length) {
				var t = isotopeSets(currTramlOut)(ichrom)
				var gc = new GhostChromatogram
				gc.precursor 		= t.q1
				gc.product 			= 0.0
				gc.collisionEnergy 	= 0.0
				gc.chromatogram.userParams += occurenceParam(t.occurence)
				gc.times 			= _times
				gc.intensities 		= chroms(ichrom)
				gc.timeDef 			= GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_LINEAR, GhostBinaryDataArray.Time(), false)
				gc.intensityDef 	= GhostBinaryDataArray.DataDef(true, false, 
											MSNumpress.ACC_NUMPRESS_SLOF, GhostBinaryDataArray.Intensity(), false)
				gc.toChromatogram(index).write(w, null)
				index += 1
			}
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
