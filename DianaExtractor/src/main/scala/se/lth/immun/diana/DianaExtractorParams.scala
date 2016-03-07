package se.lth.immun.diana

import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.BufferedReader
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream

import se.jt.Params
import se.lth.immun.xml.XmlReader

class DianaExtractorParams extends Params {

	import Params._
	
	val mzML 	= ReqString("The DIA/SWATH file to analyze (mzML)")
	val tramls 	= ReqString("The assays to extract (traml files)")
	
	val msResolution = "10"  	## "Interpreted as ppm if integer, and as m/z if float"
	val mode = "uniform"		## "Mode for XIC extraction, best|uniform|normal|square|tri"
	val diaWindows = "gillet"	## "Mode of dia window detection, gillet | snoop(guess from file) (default: gillet (1 MS1 + 32 MS2))"
	val force = false			## "Ignore missing attributes in mzML files"
	val nIsotopes = 0			## "The number highest precursor isotopes to extract"
	val noBar = false			## "Don't print process bar"
	val subSampleMs1 = 1		## "Sum k consequent ms1 spectra and store as one value"
	val outDir = "."			## "Directory where chromatogram mzML files should be stored."
	
	
	def makeBar = !noBar
	
	trait MzThreshold { def cutoff(mz:Double):Double }
	case class AbsoluteThreshold(diff:Double) extends MzThreshold { 
		override def toString = diff+" Da" 
		def cutoff(mz:Double) = diff
	}
	case class PPMThreshold(diff:Double) extends MzThreshold {
		override def toString = diff+" ppm"
		def cutoff(mz:Double) = mz / 1e6 * diff
	}
	
	lazy val mzThreshold =
		try {
			PPMThreshold(msResolution.value.toInt)
		} catch {
			case _:Throwable => 
				AbsoluteThreshold(msResolution.value.toDouble)
		}
		
	lazy val tramlFiles =
		tramls.split(" ").map(path => new File(path))
		
	lazy val mzMLFile = new File(mzML.value)
	def mzMLReader =
		new XmlReader(
			if (mzML.value.toLowerCase.endsWith(".gz"))
				new BufferedReader(new InputStreamReader(
					new GZIPInputStream(new FileInputStream(mzMLFile))))
			else
				new BufferedReader(new FileReader(mzMLFile))
		)
}