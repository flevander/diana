package se.lth.immun.diana

import se.jt.Params
import scala.util.{Try, Success, Failure}

object TracerParams {
	val UNIFORM = 0
	val BEST = 1
	val NORMAL = 2
	val SQUARE = 3
	val TRI = 4
}

class TracerParams(val name:String, val version:String) extends Params {

	import Params._
	import TracerParams._
	
	var minDiffCutoff			= 0.0
	var minDiffPPM				= 10
	var imode					= 0
	
	val mzML = ReqString("The mzML file")
	val tramls = ReqList("1-n tramls to trace")
	
	val msResolution = minDiffPPM.toString	## "Interpreted as ppm if integer, and as m/z if float"
	val mode = "uniform" ## "Mode for XIC extraction, best|uniform|normal|square|tri (default: uniform)"
	val diaWindows = "snoop" ## "Mode of dia window detection, gillet | snoop(guess from file)"
	val force = false ## "Ignore missing attributes in mzML files"
	val verbose = false ## "produce a lot more output"
	val subSampleMs1 = 1 ## "Sum k consequent ms1 spectra and store as one value"
	val outDir = "-" ## "Directory where chromatogram mzML files should be stored."
	val profiling = false ## "turn on CPU profiling"
	val timeFixedPoint = 100000.0 ## "Fixed point to used for numLin compress of time arrays"
	val intFixedPoint = 22.0 ## "Fixed point to used for numSlof compress of intensity arrays"
	
	
	def parseMinDiff = 
		Try(
			minDiffPPM = msResolution.value.toInt
		).orElse(Try({
			minDiffCutoff = msResolution.value.toDouble
			minDiffPPM = 0
		})) match {
			case Success(u) => Nil
			case Failure(msg) => List(msg.toString)
		}
	
	def parseMode:Seq[String] = {
		imode = mode.value match {
			case "uniform" => UNIFORM
			case "best" => BEST
			case "normal" => NORMAL
			case "square" => SQUARE
			case "tri" => TRI
			case _ => return List("Unknown trace mode '"+mode.value+"'")
		}
		return Nil
	}
}