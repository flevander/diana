package se.lth.immun.diana

import java.io.File
import se.jt.Params


class DianaScorerParams(val name:String, val version:String) extends Params {
	
	import Params._
	
	// USER EXPOSED PARAMS
	val verbose 		= false 	## "increase details in output"
	val concurrency 	= 1 		## "the number of assays to analyze in parallel"
	val profiling 		= false 	## "set to enable CPU profiling"
	val outEsv 			= true 		## "if false the only output is disabled"
	val output 			= "" 		## "file where result should be saved (default: input chrom .xml/.esv)"
	val pCutoff 		= 0.99 		## "p-value cutoff for peak candidates (default: 0.99)"
	val nReport			= 10		## "number of random assay to export control figure for"
	val reportSeed		= -1L		## "seed to use for report assay selection (<0 means random)"
	
	val traml = ReqString("The assay TraML file to analyze")
	val swathChromMzML = ReqString("The chromatogram MzML file to analyze")
	
	// INTERNAL PARAMS
	val relRatioUpperBound = 1.5
	val signalBinSize = DianaUtil.DEFAULT_BIN_SIZE // 20
	
	// INTERNAL VARS
	var swathFile:File = _
	var tramlFile:File = _
	var outFile:File = _
	
	var startTime		= 0L
	var chromReadTime 	= 0L
	var tramlReadTime 	= 0L
	var esvWriteTime 	= 0L
}