package se.lth.immun.diana

import scala.collection.mutable.ArrayBuffer
//import scala.collection.JavaConversions._

import swing._
import swing.event._ 

import java.awt.Color
import java.io.File
import java.io.FileReader
import java.io.BufferedReader

import se.lth.immun.math.Matrix
import se.lth.immun.math.Stats
import se.lth.immun.math.Histogram
import se.lth.immun.signal.Filter
import se.lth.immun.app.CLIApplication
import se.lth.immun.app.CommandlineArgumentException
import se.lth.immun.xml.XmlReader
import se.lth.immun.mzml.ghost._
import se.lth.immun.traml.ghost._

import se.lth.immun.math.Ratios

import se.lth.immun.anubis.ReferenceFile
import se.lth.immun.anubis.ReferencePrecursor

import se.lth.immun.anubis.Ratio
import se.lth.immun.anubis.Peak
import se.lth.immun.anubis.AnubisParams
import se.lth.immun.anubis.AnubisInput
import se.lth.immun.anubis.PepQuant
import se.lth.immun.anubis.PeakCandidate
import se.lth.immun.anubis.RelationInterval

import se.lth.immun.graphs.LineGraph
import se.lth.immun.graphs.OverlayGraph
import se.lth.immun.graphs.OnOffGraph
import se.lth.immun.graphs.OnOffCurve
import se.lth.immun.graphs.OnOffSlot
import se.lth.immun.graphs.OnOffSlotConnector
import se.lth.immun.graphs.GradientRenderer
import se.lth.immun.graphs.util._
import se.lth.immun.graphs.event.ZoomChanged
import se.lth.immun.graphs.swing.JTComboBox

import se.lth.immun.signal.Wavelet
import se.lth.immun.signal.MODWT
import se.lth.immun.signal.WaveletFilter
import se.lth.immun.markov.MarkovChainDistribution

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.stat.StatUtils

import se.lth.immun.chem._
import se.lth.immun.esv._

object DianaView extends SimpleSwingApplication with CLIApplication {

	class Isotope(
			val q1:Double,
			val occurence:Double
	) {}
	
	class ManualResult(
			val q1:Double,
			val sequence:String,
			var area:Double,
			val rt:Double
	) {}
	
	import WaveletFilter.Type._
	import DianaPeakCandidate._
	
	
	// anubis
	var nullDistTF		= new TextField("11")
	var transLimitTF	= new TextField("6")
	var pValueTolTF		= new TextField("0.01")
	var peakWidthTF		= new TextField("0.1")
	var aIn:AnubisInput = null
	var pepQuant 	= new PepQuant
	
	
	
	// files
	var mzMLFile:File 	= null
	var xmzML:XMzML		= null
	var refFile:File	= null
	var ref:ReferenceFile = null
	
	
	
	// ratio distributions
	var rDistribution:Array[Array[LineGraph]] = null
	var mDistribution:Array[Array[LineGraph]] = null
	var colDistributions	= new ArrayBuffer[LineGraph]
	var rowDistributions	= new ArrayBuffer[LineGraph]
	var labels	= new ArrayBuffer[Label]
	
	
	
	// denoiseing
	val DERIVATE 		= "derivate"
	val WAVELET_SMOOTH 	= "wavelet smooth"
	val PEAK_CAND 		= "peak candidates"
	var denoiseGraph 	= new LineGraph
	var chromCB:JTComboBox[XChromatogram] = new JTComboBox[XChromatogram]("chromatogram to view")
    var denoiseCB 		= new ComboBox[String](Array(PEAK_CAND, WAVELET_SMOOTH, DERIVATE))
	var waveletCB 		= new ComboBox[WaveletFilter.Type](
								Array(HAAR, D4, D6, D8, LA8, LA16, C6))
	var currChromatogram:XChromatogram = null
	var waveletLevelS	= new Slider
	var thresholdS		= new Slider
	var numBaseLineS	= new Slider
	var kBaseLineS		= new Slider
	var waveletB 		= Button("new wavelet"){ denoise(chromCB.selection.item.get) }
    
	
	
	// chrom and ratios
	var onOffGraph:OnOffGraph[Double]	= new OnOffGraph[Double] {
		var xAxis:Axis[Double] = new LinearAxis(0.0, 1.0) { 
			renderer = new LinearAxisRenderer(2) 
		}
		var xMin = xAxis.min
		var xMax = xAxis.max
	}
	var chromGraphOverlay		= new LineGraph
    var chromGraph 				= new LineGraph with OverlayGraph[Double, Double, Double] {
										val overlayGraph = chromGraphOverlay
									}
	var ratioCorrectnessTF				= new TextField("10")
    
	
	
	// pc analysis
    var pcOnOffGraph:OnOffGraph[Double]	= new OnOffGraph[Double] {
		var xAxis:Axis[Double] = new LinearAxis(0.0, 1.0) { 
			renderer = new LinearAxisRenderer(2) 
		}
		var xMin = xAxis.min
		var xMax = xAxis.max
	}
	var fragGraph 		= new LineGraph
	var isotopeGraph 	= new LineGraph
	var estimateGraph 	= new LineGraph
	var gradientGraph 	= new LineGraph
	var pcList			= new ListView[Carrier]
	var currPeakCandidate:Carrier = null
	//var pcFinder:DerivatePCFinder = new DerivatePCFinder(1.5, 20)

	var pValueCutoff 	= 0.99
	val signalP 		= DianaSignalProcessor.getDefault
	
	val fragmentState 		= new DianaChromatogramState(1.5)
	val fragmentEvaluator 	= new DianaPCEvaluator(fragmentState)
	var findEdges			= fragmentEvaluator.findEdges _
	
	val isotopeState 		= new DianaChromatogramState(1.5)
	val isotopeEvaluator 	= new DianaPCEvaluator(isotopeState)
	var nIsotopes	= 0
	var manualResultFile:File 		= null
	var manualResultEsv:EsvReader 	= null
	var manualResults				= new ArrayBuffer[ManualResult]
	
	
	
	// swing
    var mainFrame:MainFrame 	= null
    var precursorList			= new ListView[ReferencePrecursor]
    var statusL					= new Label
    var filterCB				= new CheckBox("filtered")
    var peaksCB					= new CheckBox("show peaks")
    var earlyExitCB				= new CheckBox("no fdr calc")
	var refGraph 				= new LineGraph
    
	
	
    precursorList.renderer 	= ListView.Renderer(pc => pc.mz.toInt + " "+pc.measuredTransitions.length+" " +pc.peptideSequence)
	pcList.renderer 		= ListView.Renderer(c => "%d-%d pn=%.1e cs=%.3f".format(c.g.istart, c.g.iend, c.fragmentMarkovPcsRatioProb, c.fragmentCorrScore))
	chromGraph.preferredSize 			= new Dimension(1000, 1000)
    chromGraph.style.annotColBackground = new Color(0.5f, 0.5f, 0.5f, 0.5f)
    chromGraph.cacheRender 				= true
    //chromGraphOverlay.cacheRender 	= true
    
    chromCB.minimumSize = new Dimension(300, 20)
    chromCB.preferredSize = new Dimension(300, 20)
    
    onOffGraph.connectorVisZoomFactor	= 0.01
	onOffGraph.preferredSize 			= new Dimension(1000, 1000)
    onOffGraph.cacheRender 				= true
    onOffGraph.style.backgroundColor 	= Color.LIGHT_GRAY
    
    denoiseGraph.cacheRender			= true
    denoiseGraph.style.backgroundColor 	= Color.LIGHT_GRAY
    
    pcOnOffGraph.connectorVisZoomFactor	= 0.01
	pcOnOffGraph.preferredSize 			= new Dimension(1000, 1000)
    pcOnOffGraph.cacheRender 			= true
    pcOnOffGraph.style.backgroundColor 	= Color.LIGHT_GRAY
    
    fragGraph.preferredSize 			= new Dimension(1000, 1000)
    fragGraph.style.annotColBackground 	= new Color(0.5f, 0.5f, 0.5f, 0.5f)
    fragGraph.cacheRender 				= true
    
    isotopeGraph.preferredSize 				= new Dimension(1000, 1000)
    isotopeGraph.style.annotColBackground 	= new Color(0.5f, 0.5f, 0.5f, 0.5f)
    isotopeGraph.cacheRender 				= true
    
    estimateGraph.preferredSize 			= new Dimension(200, 200)
    estimateGraph.style.annotColBackground 	= new Color(0.5f, 0.5f, 0.5f, 0.5f)
    estimateGraph.cacheRender 				= true
    
    gradientGraph.preferredSize 			= new Dimension(1000, 1000)
    gradientGraph.style.annotColBackground 	= new Color(0.5f, 0.5f, 0.5f, 0.5f)
    gradientGraph.cacheRender 				= true
    gradientGraph.renderer 					= new GradientRenderer[Double]
    
    refGraph.style.annotColBackground 	= new Color(0.5f, 0.5f, 0.5f, 0.5f)
    refGraph.preferredSize 				= new Dimension(200, 200)
    
    waveletLevelS.orientation = Orientation.Horizontal
    waveletLevelS.max = 6
    waveletLevelS.min = 1
    waveletLevelS.majorTickSpacing = 1
    waveletLevelS.paintLabels = true
    waveletLevelS.paintTicks = true
    
    thresholdS.orientation = Orientation.Horizontal
    thresholdS.max = 100
    thresholdS.min = 0
    thresholdS.value = 10
    thresholdS.majorTickSpacing = 20
    thresholdS.paintLabels = true
    thresholdS.paintTicks = true
    
    numBaseLineS.orientation = Orientation.Horizontal
    numBaseLineS.max = 100
    numBaseLineS.min = 1
    numBaseLineS.value = 40
    numBaseLineS.majorTickSpacing = 20
    numBaseLineS.paintLabels = true
    numBaseLineS.paintTicks = true
    
    kBaseLineS.orientation = Orientation.Horizontal
    kBaseLineS.max = 100
    kBaseLineS.min = 0
    kBaseLineS.value = 50
    kBaseLineS.majorTickSpacing = 20
    kBaseLineS.paintLabels = true
    kBaseLineS.paintTicks = true
    
    earlyExitCB.selected = true
	
	listenTo(chromGraph)
	listenTo(fragGraph)
	listenTo(precursorList.selection)
	listenTo(pcList.selection)
	listenTo(chromCB.selection)
	//listenTo(waveletCB.selection)
	listenTo(denoiseCB.selection)
	reactions += {
		case zc:ZoomChanged[Double] => {
			if (zc.source == chromGraph) {
				onOffGraph.setZoom(zc.start, zc.end)
				onOffGraph.repaint
			}
			if (zc.source == fragGraph) {
				pcOnOffGraph.setZoom(zc.start, zc.end)
				pcOnOffGraph.repaint
				isotopeGraph.setZoom(zc.start, zc.end)
				isotopeGraph.repaint
			}
		}
		case sc:SelectionChanged => {
			if (sc.source == precursorList) 
				if (!precursorList.selection.items.isEmpty)
					analyzeChrom(precursorList.selection.items(0))
			if (sc.source == pcList)
				if (!pcList.selection.items.isEmpty)
					highLightPCGroup(pcList.selection.items(0))
			if (sc.source == chromCB || sc.source == denoiseCB || sc.source == waveletCB)
				chromCB.selection.item match {
					case Some(xc) => denoise(xc)
					case None => clearDenoise
				}
		}
	}
    
	
    class ColorLabel(col:Color) extends Label {
    	background = col
    }
	
	
	
	override def main(args:Array[String]):Unit = {
		
		arg("REFERENCE_FILE", s => {
	    		refFile = new File(s)
    			var r = new XmlReader(new BufferedReader(new FileReader(refFile)))
    			if (s.toLowerCase.endsWith(".traml"))
    				ref = ReferenceFile.fromTraML(r)
	    		else 
	    			ref = ReferenceFile.fromFile(r)
	    	})
	    
		arg("MZML_FILE", s => {
			mzMLFile = new File(s)
			xmzML = XMzML.fromFile(new XmlReader(new BufferedReader(new FileReader(mzMLFile))))
		})
			
		opt("isotopes", 
				"The number highest precursor isotopes to extract (default 0)", 
				s => nIsotopes = s.toInt, 
				"X"
			)
			
		opt("manual", 
				"The number highest precursor isotopes to extract (default 0)", 
				s => {
					manualResultFile 	= new File(s)
					manualResultEsv		= new EsvReader(new BufferedReader(new FileReader(manualResultFile)))
					while (!manualResultEsv.EOF) {
						val q1 = manualResultEsv.getValue("PrecursorMz").toDouble
						val seq = manualResultEsv.getValue("PeptideSequence")
						manualResults.find(mr => mr.q1 == q1 && mr.sequence == seq) match {
							case Some(mr) => mr.area += manualResultEsv.getValue("Area").toDouble
							case None => 
								manualResults += new ManualResult(q1, seq,
									manualResultEsv.getValue("Area").toDouble,
									manualResultEsv.getValue("RetentionTime").toDouble
								)
						}
						manualResultEsv.readLine
					}
				}, 
				"ESV"
			)
		
		try {
			parseArgs("DianaView", args)
		} catch {
			case e:CommandlineArgumentException => 
				return
		}
		
		precursorList.listData = ref.precursors.sortBy(_.mz)
		rDistribution = Matrix.get2d[LineGraph](10)
		for (i <- 0 until 10)
			for (j <- (i+1) until 10)
				rDistribution(i)(j) = new LineGraph

		mDistribution = Matrix.get2d[LineGraph](10)
		for (i <- 0 until 10)
			for (j <- (i+1) until 10)
				mDistribution(i)(j) = new LineGraph {
					xAxis = new LogAxis(java.lang.Double.MIN_NORMAL, 1.0) { renderer = new LinearAxisRenderer(2) }
					xMin = xAxis.min
					xMax = xAxis.max
				}
		
		
		super.main(args)
	}
	
	
	
	
	def top = {
		mainFrame = new MainFrame {
			title = "DianaView - " + mzMLFile
			
			contents = new BorderPanel {
				focusable = true
				minimumSize = new Dimension(1024, 800)
				preferredSize = new Dimension(1024, 800)
				
				import BorderPanel.Position._
				
				layout(new BoxPanel(Orientation.Vertical) {
					contents += new ScrollPane(precursorList)
					contents += new ScrollPane(pcList)
					contents += refGraph
					contents += estimateGraph
				}) = West
				layout(new BorderPanel {
					layout(new GridPanel(2, 6) {
						contents += new Label("null dist")
						contents += nullDistTF
						contents += new Label("trans limit")
						contents += transLimitTF
						contents += new Label("...")
						contents += peaksCB
						contents += new Label("peak width")
						contents += peakWidthTF
						contents += filterCB
						contents += earlyExitCB
						contents += new Label("ratio correctness scaling")
						contents += ratioCorrectnessTF
					}) = North
					layout(new TabbedPane() {
						pages += new TabbedPane.Page(
								"pc analysis",
								new BoxPanel(Orientation.Vertical) {
									contents += fragGraph
									contents += isotopeGraph
									contents += pcOnOffGraph//gradientGraph
								})
						pages += new TabbedPane.Page(
								"chromatogram and ratios", 
								new BoxPanel(Orientation.Vertical) {
									contents += chromGraph
									contents += onOffGraph
								})
						pages += new TabbedPane.Page(
								"distributions", 
								new GridPanel(11, 11) {
									var l:LineGraph = null
									var label:Label = null
									for (r <- 0 until 11)
										for (c <- 0 until 11) {
											if (c == 0 && r == 0) //TOP LEFT
												contents += new Label
											else if (c > 0 && c == r) {//DIAGONAL
												label = new Label(""+c)
												labels		+= label
												contents 	+= label
											} else if (r == 0 && c > 0) {//COL HEADERS
												l = new LineGraph
												l.style.backgroundColor = Color.LIGHT_GRAY
												l.style.curveColors = Array(Color.CYAN)
												colDistributions += l
												contents += l
											} else if (c == 0 && r > 0) {//ROW HEADERS
												l = new LineGraph
												l.style.backgroundColor = Color.LIGHT_GRAY
												l.style.curveColors = Array(Color.CYAN)
												rowDistributions += l
												contents += l
											} else if (c > r) //DISTRIBUTIONS
												contents += rDistribution(r-1)(c-1)
											else //REST
												contents += mDistribution(c-1)(r-1)
										}
								})
						pages += new TabbedPane.Page(
								"denoising",
								new BorderPanel {
									layout(new FlowPanel {
										contents += chromCB
										contents += denoiseCB
										contents += waveletCB
										contents += thresholdS
										contents += waveletLevelS
										contents += numBaseLineS
										contents += kBaseLineS
										contents += waveletB
									}) = North
									layout(denoiseGraph) = Center
								})
					}) = Center
				}) = Center
				layout(statusL) = South
			}
		}
		mainFrame
	} // top
	
	
	
	
	
	var currPrecursor:Double = -1.0
	def analyzeChrom(pc:ReferencePrecursor):Unit = {
		if (pc.mz == currPrecursor) return
		pepQuant.earlyExit = earlyExitCB.selected
		var params = new AnubisParams(
						peakWidthTF.text.toDouble,
						nullDistTF.text.toInt,
						transLimitTF.text.toInt,
						false,
						pValueTolTF.text.toDouble)
		
		try {
        	xmzML.grouper.extractGroup(pc.mz) match {
            	case Some(cg1) => {
            		aIn = AnubisInput.of(cg1, pc, params)
            		var peaks = pepQuant.findBestCandidate(aIn)
            		
            		println("########## " +pc.peptideSequence+ " ###########")
            		println("   ratio calc: "+pepQuant.t0 + "ms")
            		println("     find pcs: "+pepQuant.t1 + "ms")
            		println(" sanity check: "+pepQuant.t2 + "ms")
            		println("null dist gen: "+pepQuant.t3 + "ms")
            		println("     s/n calc: "+pepQuant.t4 + "ms")
            		println("     fdr calc: "+pepQuant.t5 + "ms")
            		println("        quant: "+pepQuant.t6 + "ms")
            		
            		var now = System.currentTimeMillis
            		//setChrom(pc.mz)
            		var iAIn:AnubisInput = null
            		if (nIsotopes > 1) {
            			iAIn = getIsotopeInput(xmzML, pc)
            			showIsotopes(iAIn.cg)
            		}
            		showRef(aIn)
            		//display(aIn.cg, peaks)
            		//showRatios(aIn.cg)
            		//showDistributions
            		chromCB.model.setItems(aIn.cg.chromatograms)
            		pcAnalysis(aIn, iAIn, pc)
            		currPrecursor = pc.mz
            		println("      display: "+(System.currentTimeMillis - now) + "ms")
            	}
            	case None => {
            		status(pc.proteinName+" "+pc.peptideSequence+"{M/Z: "+pc.mz+"} not found in file '"+mzMLFile+"'")
            	}
            }
        }
        catch {
        	case e:Exception => {
        		status("YARGH! "+pc.proteinName+" "+pc.peptideSequence+"{M/Z: "+pc.mz+"}")
        		e.printStackTrace
        	}
        }
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
		
		return new AnubisInput(seq, q1, cg, isotopes.map(i => AnubisInput.Fragment(i.q1, "isotope")), r.toArray, new AnubisParams)
	}
	
	
	
	
	
	def showIsotopes(cg:XChromatogramGroup) = {
		val chroms 	= cg.chromatograms 
		val l 		= chroms.length
		val t 		= chroms(0).times
		
		isotopeGraph.setCurves(chroms.reverse.map(chrom => {
				var y 	= chrom.intensities
				new Curve2(t, y, y.map(d => false), "isotope %.2f".format(chrom.q1))
			}))
		isotopeGraph.repaint
	}
	
	
	
	
	
	var currFiltering:Boolean = false
	/*
	def setChrom(q1:Double):Unit = {
		
		import anubis.Ratio
		
		if (currPrecursor == q1 && currFiltering == filterCB.selected) return
		
		if (aIn.cg.chromatograms.isEmpty) {
			chromGraph.clear
			return status("no chromatograms for precursor "+q1+"m/z")
		}
		
		var ratios = pepQuant.ratios
		state.setChromatogram(
				aIn.cg.chromatograms.length, 
				ratios, 
				ratios.getTargetRatioTable(Ratio.ONE)(aIn.refRatios, Ratio.INVERSE)
			)
		
		var times = aIn.cg.chromatograms(0).times
		var ratioCorrectnessScale = ratioCorrectnessTF.text.toDouble
		var ratioCorrectnessCurve = aIn.refRatios.map(r => {
			ratios.getRatio(r.transitionId1, r.transitionId2).map(d => {
				if (!java.lang.Double.isNaN(d) && d > 0 && d < Double.PositiveInfinity)
					1 / (1 + ratioCorrectnessScale * math.abs(math.log(d) - math.log(r.mean)))
				else
					0.0
			})
		}).transpose.map(_.sum)
		
		/*
		chromGraphOverlay.setCurves(Array(
					new Curve2(
						times, 
						jt.signal.Filter.savitzkyGolay9(ratioCorrectnessCurve), 
						ratioCorrectnessCurve.map(d => java.lang.Double.isNaN(d)),
						"ratio correctness",
						new Color(0.0f, 0.0f, 0.0f, 0.5f)
					)))
		*/
		if (filterCB.selected)
			chromGraph.setCurves(
				pepQuant.filtered.zip(aIn.cg.chromatograms).map(
					t => new Curve2(
							times,
							t._1,
							t._1.map(d => java.lang.Double.isNaN(d)),
							q1+" "+t._2.q3
						)
				).toSeq
			)
		else
			chromGraph.setCurves(
				aIn.cg.chromatograms.map(
					c => new Curve2(
							times,
							c.intensities,
							c.intensities.map(d => java.lang.Double.isNaN(d)),
							c.q1+" "+c.q3
						)
				).toSeq
			)
		
		chromGraph.title = "#datapoints: "+aIn.cg.chromatograms(0).length
		
		currFiltering = filterCB.selected
		currPrecursor = q1
		chromGraph.repaint
		
        status("loaded chromatogram "+aIn.precursorMZ+"m/z    "+pepQuant.uncheckedPcs.length+" pcs")
	}
	
	*/
	
	
	
	def display(cg:XChromatogramGroup, peaks:Seq[Peak]) = {
		val times 	= cg.chromatograms(0).times
		//pcList.listData = pepQuant.uncheckedPcs
		for (pc <- pepQuant.uncheckedPcs) {
			chromGraph.addAnnotation(
					new HeightBoxAnnotation(
							math.min(times.last, times(math.max(pc.start, 0))), 
							pc.estimates.map(_.max).max, 
							math.max(times.head, times(math.min(pc.end, times.length - 1))), 
							Annotation.BACKGROUND,
							if (pepQuant.sanePcs.contains(pc)) "fdr: "+pc.fdr
							else	"insane"
					))
		}
	}
	
	
	/*
	def highLightPC(pc:PeakCandidate):Unit = {
		if (currPeakCandidate == pc) return
		currPeakCandidate = pc
		//println
		//pc.conditions.foreach(c => println(c.ratio.ratio+" start:" +c.start+" end:"+c.end))
		
		chromGraph.annotations = chromGraph.annotations.filter(!_.active)
		onOffGraph.connectors = onOffGraph.connectors.filter(!_.active)
		
		val times 	= aIn.cg.chromatograms(0).times
		chromGraph.addAnnotation(
					new HeightBoxAnnotation(
							math.min(times.last, times(math.max(pc.start, 0))), 
							pc.estimates.map(_.max).max, 
							math.max(times.head, times(math.min(pc.end, times.length - 1))), 
							Annotation.ACTIVE,
							if (pepQuant.sanePcs.contains(pc)) "fdr: "+pc.fdr
							else	"insane"
					))
		chromGraph.repaint
					
		var sorted = pc.conditions.sortBy(_.start)
		var edges = pepQuant.pcf.findEdges(sorted)
		for (i <- 0 until edges.length)
			for (j <- edges(i)) {
				var a = sorted(i)
				var b = sorted(j)
				var x0 = math.max(times(a.start), times(b.start))
				var xn = math.min(times(a.end), times(b.end))
				var mid = (xn + x0) / 2
				//println("edge:"+mid+"   "+a.ratio.ratio+"-"+b.ratio.ratio)
				onOffGraph.connectors += new OnOffSlotConnector(a.ratio.ratio, mid, b.ratio.ratio,	mid, true)
			}
		
		onOffGraph.repaint
	}
	*/
	
	
	
	/*
	def showRatios(cg:XChromatogramGroup) = {
		val l 		= cg.chromatograms.length
		val times 	= cg.chromatograms(0).times
		var connectors = new ArrayBuffer[OnOffSlotConnector[Double]]
		val curves 	= new Array[OnOffCurve[Double]]((l*(l-1)/2))
		Ratios.iterate(l, (ri, ui, di) => {
			curves(ri) = new OnOffCurve[Double](
					new ArrayBuffer[OnOffSlot[Double]],
					cg.chromatograms(ui).q3.toInt + " / " +cg.chromatograms(di).q3.toInt,
					ui,
					di
			)
		})
		
		var i = 0
		for (pc <- pepQuant.uncheckedPcs) {
			for (ri <- pc.conditions)
				curves(ri.ratio.ratio).slots += 
					new OnOffSlot(times(ri.start), times(ri.end), None, i)
			
			i += 1
			
		}
		
		curves.foreach(c => c.slots = c.slots.sortBy(_.start))
		onOffGraph.setCurves(curves, chromGraph.xMin, chromGraph.xMax)
		onOffGraph.connectors = connectors
		onOffGraph.setZoom(chromGraph.xAxis.min, chromGraph.xAxis.max)
		onOffGraph.repaint
	}
	*/
	
	
	
	def showDistributions = {
		var numFragments = pepQuant.filtered.length
		for (i <- 0 until numFragments) {
			var xc 		= aIn.cg.chromatograms(i)
			var hist 	= Histogram.fromArray(xc.intensities.map(d => math.log10(d + 1)).toArray, 20)
			labels(i).text = "%.2f".format(xc.q3)
			colDistributions(i).setCurves(Array(
				new Curve2(hist.bins, hist.counts, hist.counts.map(_ => false), labels(i).text)
			))
			rowDistributions(i).setCurves(Array(
				new Curve2(hist.bins, hist.counts, hist.counts.map(_ => false), labels(i).text)
			))
		}
		
		var mults = new ArrayBuffer[(Array[Double], String)]
			
		Ratios.iterate(numFragments, (ri, ui, di) => {
			var r = pepQuant.ratios.getRatio(ui, di)
						.filter(d => d > 0.0 && d < Double.PositiveInfinity)
						.map(d => math.log10(d))
			var median 	= StatUtils.percentile(r, 50)
			var mad 	= Stats.mad(r)
			r = r.filter(d => math.abs(d - median) < 5*mad)
			if (!r.isEmpty) {
				var hist = Histogram.fromArray(r, 20)
				rDistribution(ui)(di).setCurves(List(
					new Curve2(hist.bins, hist.counts, hist.bins.map(_ => false), "")
				))
				aIn.refRatios.find(r => 
					r.transitionId1 == ui && r.transitionId2 == di
				) match {
					case Some(ratio) => { 
						rDistribution(ui)(di).addAnnotation(
								new XAnnotation(
										math.log10(ratio.mean * 2),
										Annotation.ACTIVE
								))
						rDistribution(ui)(di).addAnnotation(
								new XAnnotation(
										math.log10(ratio.mean / 2),
										Annotation.ACTIVE
								))
						var zero = new XAnnotation(	0.0)
						zero.background = true
						rDistribution(ui)(di).addAnnotation(zero)
					}
					case None => {}
				}
			}
			
			var y1 = aIn.cg.chromatograms(ui).intensities.toArray
			var y2 = aIn.cg.chromatograms(di).intensities.toArray
			var L	= y1.length
			var m = new Array[Double](L)
			for (i <- 0 until L)
				m(i) = y1(i) * y2(i)
			mults += new Tuple2(m, aIn.cg.chromatograms(ui).q3 +"*"+ aIn.cg.chromatograms(di).q3)
			//println("num < 0: "+m.filter(d => d < 0.0).length)
			//println("num = 0: "+m.filter(d => d == 0.0).length)
			var y1m 	= StatUtils.percentile(y1.filter(_ > 0.0), 50)
			var y2m 	= StatUtils.percentile(y2.filter(_ > 0.0), 50)
			var mid		= y1m * y2m
			//mad 	= JTMath.stdDevMAD(m)
			//m = m.filter(d => math.abs(d - median) < 5*mad)
			var hist = Histogram.fromArray(m, 20, true)
			mDistribution(ui)(di).setCurves(List(
					new Curve2(hist.bins, hist.counts, hist.bins.map(_ => false), "")
				))
			mDistribution(ui)(di).addAnnotation(
					new XAnnotation(y1m))
			mDistribution(ui)(di).addAnnotation(
					new XAnnotation(y2m))
			mDistribution(ui)(di).addAnnotation(
					new XAnnotation(mid, Annotation.ACTIVE))
		})
		
		var max = Double.MinValue
		var min = Double.MaxValue
		Ratios.iterate(numFragments, (ri, ui, di) => {
			var md = mDistribution(ui)(di) 
			if (md.xAxis.max > max) max = md.xAxis.max
			if (md.xAxis.min < min) min = md.xAxis.min  
		})
		Ratios.iterate(numFragments, (ri, ui, di) => {
			mDistribution(ui)(di).setZoom(min, max) 
		})
	}
	
	
	
	
	def denoise(xc:XChromatogram) = {
		val missing 	= xc.intensities.map(java.lang.Double.isNaN _)
		val t 			= xc.times
		val intensities = xc.intensities.toArray
		var str = denoiseCB.selection.item 
		if (str == DERIVATE) {
			var y 		= se.lth.immun.signal.Filter.savitzkyGolay9(intensities)
			var dy 		= y.zip(y.tail).map(	tu => 	 tu._2 - tu._1		)
			var dyt 	= t.zip(t.tail).map(	tu => 	(tu._2 + tu._1)/2	)
			var ddy		= dy.zip(dy.tail).map(	tu => 	 tu._2 - tu._1		)
			var ddyt 	= dyt.zip(dyt.tail).map(tu => 	(tu._2 + tu._1)/2	)
			
			denoiseGraph.setCurves(Array(
				new Curve2(t, 		intensities, missing, "original", Color.GRAY),
				new Curve2(t, 		y, 			missing, "y", Color.BLUE),
				new Curve2(dyt, 	dy, 		missing, "dy", Color.RED),
				new Curve2(ddyt, 	ddy, 		missing, "ddy", Color.ORANGE)
			))
			denoiseGraph.addAnnotation(new YAnnotation(0.0))
			status("loaded derivates "+xc.q1 +" "+xc.q3)
			
			
			
		} else if (str == WAVELET_SMOOTH) {
			var waveletLevels = waveletLevelS.value
			var filter = new WaveletFilter(waveletCB.selection.item)
			var modwt = new MODWT
			var boundary = Wavelet.Boundary.PERIODIC;
			var wavelets = Wavelet.decompose(intensities, xc.length, waveletLevels, 
												filter, modwt, boundary, null)
			
			var smoothLevels = Wavelet.smooths(wavelets, waveletLevels, filter, modwt, boundary, null);
			var curves = new Array[Curve2[Double, Double]](smoothLevels.length)
			for (i <- 0 until smoothLevels.length) {
				curves(i) = new Curve2(t, smoothLevels(i), missing, "wavelet smooth "+(1+i), new Color(0x220000 * i))
			}
			var kThreshold = (thresholdS.value) / 100.0
			var threshold = math.pow(10, kThreshold * 3.0 - 2.0)
			var denoised = Wavelet.denoise(wavelets, waveletLevels, filter, modwt, threshold, null)
			
			var noBaseLine 	= numBaseLineS.value
			var kBaseLine 	= kBaseLineS.value * 0.01
			var binSize 	= t.length / noBaseLine
			var baseline = t.zip(
								smoothLevels(waveletLevels-1).map(d => math.max(0.0, d))
							).grouped(binSize).map(g => 
								g.min(Ordering[Double].on[(Double, Double)](_._2))
							).toSeq
			var baseline2 = xc.intensities.sliding(binSize).map(w => StatUtils.percentile(w.toArray, 50)).toSeq
			var baseline3 = smoothLevels(waveletLevels-1).sliding(binSize).map(w =>{ 
								val sort = w.sorted
								val bottom = sort.take((w.length * kBaseLine).toInt)
								val mad = Stats.mad(bottom)
								StatUtils.mean(bottom) + kThreshold*mad
							}).toSeq
							
			denoiseGraph.setCurves(
					Array(
							new Curve2(t, xc.intensities, missing, "original", Color.BLUE),
							new Curve2(t, denoised, missing, "denoised: "+threshold, Color.YELLOW),
							curves(waveletLevels-1),
							new Curve2(baseline.map(_._1), baseline.map(_._2), baseline.map(_ => false), "baseline", Color.WHITE),
							new Curve2(t.drop(binSize / 2).take(baseline2.length), baseline2, baseline2.map(_ => false), "baseline2", Color.CYAN),
							new Curve2(t.drop(binSize / 2).take(baseline3.length), baseline3, baseline3.map(_ => false), "baseline3", Color.PINK)
					))
			
			denoiseGraph.addAnnotation(new YAnnotation(0.0) {background = true})
			/*									
			var curves = Array( 
					new Curve2(t, wsmoothed.last.smooth, missing, "wavelet "+waveletLevels, new Color(0x220000 * waveletLevels)),
					new Curve2(t, xc.intensities, missing, "original data", Color.BLUE)
				)
			denoiseGraph.setCurves(curves)
			*/
			/*
			var wsmoothed = Wavelet.decompose(xc.intensities, xc.length, waveletLevels, 
												filter, modwt, boundary, null, true)
			
												
			println(xc.intensities.take(12).mkString(" "))
			println(missing.take(12).mkString(" "))
			wsmoothed.foreach(w => println(w.smooth.take(12).mkString(" ")))
			var curves = wsmoothed.map(w => 
					new Curve2(t, w.smooth, missing, "wavelet "+w.level, new Color(0x220000 * w.level))
				)
			denoiseGraph.setCurves(
				new Curve2(t, xc.intensities, missing, "original", Color.BLUE) +: curves						
			)
			*/
			status("loaded wavelet "+xc.q1 +" "+xc.q3+" with "+waveletCB.selection.item+" wavelet")
		} else if (str == PEAK_CAND) {
			var filter = new WaveletFilter(waveletCB.selection.item)
			var modwt = new MODWT
			var boundary = Wavelet.Boundary.PERIODIC;
			var wavelets = Wavelet.decompose(xc.intensities.toArray, xc.length, 2, 
												filter, modwt, boundary, null)
			var smoothLevels = Wavelet.smooths(wavelets, 2, filter, modwt, boundary, null)
			var y		= smoothLevels(1)
			var dy 		= y.zip(y.tail).map(	tu => 	 tu._2 - tu._1		)
			var dyt 	= t.zip(t.tail).map(	tu => 	(tu._2 + tu._1)/2	)
			var ddy		= dy.zip(dy.tail).map(	tu => 	 tu._2 - tu._1		)
			var ddyt 	= dyt.zip(dyt.tail).map(tu => 	(tu._2 + tu._1)/2	)
			
			var noBaseLine 	= numBaseLineS.value
			var binSize 	= t.length / noBaseLine
			var baseline2 = y.sliding(binSize).map(w => math.max(1.0, StatUtils.percentile(w, 50))).toSeq
			
			denoiseGraph.setCurves(
					Array(
							new Curve2(t, xc.intensities, missing, "original", Color.GRAY),
							new Curve2(t.drop(binSize / 2).take(baseline2.length), baseline2, baseline2.map(_ => false), "baseline", Color.WHITE),
							new Curve2(t, y, missing, "y: wavelet smooth 2", Color.BLUE),
							new Curve2(dyt, dy, missing, "dy", Color.CYAN),
							new Curve2(ddyt, ddy, missing, "ddy", Color.PINK)
					))
					
			var pcs = DianaPeakCandidate.findPCs(y, dy, ddy, baseline2, binSize)
			for (pc <- pcs)
				denoiseGraph.addAnnotation(
						new HeightBoxAnnotation(
							t(pc.istart), 
							y(pc.iapex), 
							t(pc.iend), 
							Annotation.BACKGROUND
					))

			denoiseGraph.addAnnotation(new YAnnotation(0.0) {background = true})
			status("loaded peak candidates for "+xc.q1 +" "+xc.q3)
		}
		
		denoiseGraph.title = xc.q1 +" "+xc.q3
		denoiseGraph.repaint
	}
	
	
	
	
	def pcAnalysis(aIn:AnubisInput, iAIn:AnubisInput, pc:ReferencePrecursor):Unit = {
		
		val before 		= System.currentTimeMillis
		
		val fChroms = aIn.cg.chromatograms
		val iChroms:Seq[XChromatogram] = if (iAIn != null) iAIn.cg.chromatograms else Nil
		
		val allChroms 	= fChroms ++ iChroms
		if (allChroms.isEmpty) return
		val l 		= fChroms.length
		val t 		= fChroms(0).times
		val curves 	= new Array[OnOffCurve[Double]](l)
		
		for (i <- 0 until l)
			curves(i) = new OnOffCurve[Double](
					new ArrayBuffer[OnOffSlot[Double]],
					fChroms(i).q1.toInt + " / " +fChroms(i).q3.toInt,
					i, i
				)
		
		
		def getPCs(y:DianaSignalProcessor.SmoothAndBase) = {	
			var dy 		= signalP.getDerivate(y.smooth)
			var ddy		= signalP.getDerivate(dy)
			
			DianaPeakCandidate.findPCs(y.smooth, dy, ddy, y.base)
		}
		
		var smooths		= allChroms.map(xc => signalP.getSmoothAndBase(xc.intensities.toArray))
		var pcs 		= smooths.map(xc => getPCs(xc))
		var carriers 	= DianaPeakCandidate.groupPCs(pcs, findEdges)
							.map(g => new Carrier(pc, g, fChroms.length, iChroms.length))
							.filter(_.fragmentPcs.length > 0)
		
		var reduced = Filter.baseLineReduce(
        					fChroms.toArray.map(tr => tr.intensities.toArray)
        				)
        var filtered 	= reduced.map(ds => Filter.savitzkyGolay9(ds))
        var ratios 		= new Ratios(filtered)
		fragmentState.setChromatogram(
				l, 
				ratios, 
				ratios.getTargetRatioTable(Ratio.ONE)(aIn.refRatios, Ratio.INVERSE),
				pc.retentionTime.peak,
				t.toArray
			)
			
		var i = 0
		for (c <- carriers) {
			if (c.fragmentPcs.length > 1) {
				for (pc <- c.fragmentPcs)
					curves(pc.icurve).slots += new OnOffSlot(t(pc.istart), t(pc.iend), Some(t(pc.iapex)), i)
				i += 1
				c.fragmentValids = fragmentEvaluator.validateGroup(c.g, c.fragmentPcs, smooths)
			} else {
				var pc = c.fragmentPcs.head
				curves(pc.icurve).slots += new OnOffSlot(t(pc.istart), t(pc.iend), Some(t(pc.iapex)))
			}
		}
		
		carriers = carriers.filter(_.fragmentValids != null)
		if (carriers.nonEmpty)
			fragmentState.calculateChromatogramStatistics(carriers.map(c => (c.g, c.fragmentValids)))
		
		for (c <- carriers) {
			c.fragmentRankAllRatioProb = fragmentEvaluator.nullRatioProb(
									c.g, 
									c.fragmentValids,
									fragmentState.statsAll,
									"rank")
			c.fragmentRankPcsRatioProb = fragmentEvaluator.nullRatioProb(
									c.g, 
									c.fragmentValids,
									fragmentState.statsPcs,
									"rank")
			c.fragmentMarkovAllRatioProb = fragmentEvaluator.nullRatioProb(
									c.g, 
									c.fragmentValids,
									fragmentState.statsAll,
									"markov")
			c.fragmentMarkovPcsRatioProb = fragmentEvaluator.nullRatioProb(
									c.g, 
									c.fragmentValids,
									fragmentState.statsPcs,
									"markov")
		}
		
		
		pcOnOffGraph.setCurves(curves, t.head, t.last)
		pcOnOffGraph.setZoom(t.head, t.last)
		pcOnOffGraph.repaint
		
		pcList.listData = carriers.sortBy(_.fragmentMarkovPcsRatioProb)
		
		fragGraph.setCurves(smooths.zip(fChroms).map(yxc => {
				var y 	= yxc._1.smooth
				var xc 	= yxc._2
				new Curve2(xc.times, y, y.map(d => false), "fragment "+xc.q3.toInt)
			}))
		if (peaksCB.selected)
			for (c <- carriers) {
				fragGraph.addAnnotation(
						new HeightBoxAnnotation(
								t(math.max(0, math.min(c.g.istart, t.length - 1))), 
								c.fragmentPcs.map(pc => smooths(pc.icurve).smooth.slice(pc.istart, pc.iend).max).max, 
								t(math.max(0, math.min(c.g.iend, t.length - 1))), 
								Annotation.BACKGROUND,
								"p:%.1e".format(c.fragmentMarkovPcsRatioProb)
						))
			}
	
		/*
		if (!validGroups.isEmpty) {
			var best = validGroups.min(Ordering[Double].on[PCGroup](_.pvalue))
			println("         best: "+best+"    p="+best.pvalue+"    area="+best.area+"   pEst="+best.pEstimates+"    corrScore="+best.corrScore)
			status("chromatogram for "+aIn.peptideSequence+" analyzed, best candidate found with p="+best.pvalue+" and area="+best.area)
			
			for (g <- validGroups) {
				pcChromGraph.addAnnotation(
						new HeightBoxAnnotation(
								t(math.max(0, math.min(g.istart, t.length - 1))), 
								g.pcs.map(pc => smooths(pc.icurve).smooth.slice(pc.istart, pc.iend).max).max, 
								t(math.max(0, math.min(g.iend, t.length - 1))), 
								if (g.pvalue == best.pvalue && g.pvalue < 0.05) Color.RED else Annotation.BACKGROUND,
								"p:%.1e".format(g.pvalue)
						))
			}
		
			showEstimates(best)
		}
		*/
		if (manualResults.nonEmpty) {
			manualResults.find(mr => 
						math.abs(mr.q1 - pc.mz) < 0.01
					&& 	mr.sequence == pc.peptideSequence.takeWhile(_ != '(')
			) match {
				case Some(mr) =>
					fragGraph.addAnnotation(new XAnnotation(mr.rt, Color.CYAN, "manual"))
				case None => {}
			}
		}
		fragGraph.addAnnotation(new XAnnotation(pc.retentionTime.peak * 0.55741152 + 36.10350925, Color.BLUE, "lib rt"))
		fragGraph.setZoom(t.head, t.last)
		fragGraph.repaint
	}
	
	
	
	
	def highLightPCGroup(c:Carrier):Unit = {
		if (currPeakCandidate == c) return
		currPeakCandidate = c
		//println
		//pc.conditions.foreach(c => println(c.ratio.ratio+" start:" +c.start+" end:"+c.end))
		
		fragGraph.annotations = fragGraph.annotations.filter(!_.active)
		
		val times 	= aIn.cg.chromatograms(0).times
		for (a <- fragGraph.annotations) {
			a match {
				case h:HeightBoxAnnotation[Double, Double] => 
					if (h.x1 == times(math.min(times.length-1, math.max(c.g.istart, 0)))
						&& 	h.x2 == times(math.min(times.length-1, math.max(c.g.iend, 0)))
					) fragGraph.addAnnotation(new HeightBoxAnnotation(
							h.x1, h.y1, h.x2, Annotation.ACTIVE, "p:%.1e".format(c.fragmentMarkovPcsRatioProb)))
				case _ => {}
			}
		}
		
		//fragmentEvaluator.pvalueForFragments(g, aIn.cg.chromatograms.map(_.intensities))
		/*
		println("PC     start:%.1f    width:%d     total p: %.1e    pEst: %.1e   corr: %.3f".format(
				times(math.min(times.length-1, math.max(g.istart, 0))), 
				g.iend - g.istart,
				g.pvalue,
				g.pEstimates,
				g.corrScore
			))
		println("fragment p values: "+g.fragmentPvalues.map(p => "%.1e".format(p)).mkString(" "))
		println("      fragment xs: "+g.xs.map(p => "%.1f".format(p)).mkString(" "))
		println("  fragment sigmas: "+g.sigmas.map(p => "%.1f".format(p)).mkString(" "))
		*/
			/*
		println(g.pvalueTable.map(
					_.map(p => "%.2e".format(p)).mkString(" ")
				).mkString(" "))
				*/
		
		println("VALIDATIONS: ")
		c.fragmentValids.valids.foreach(v => println(aIn.cg.chromatograms(v.icurve).q3.toInt + 
											"       \t" +
											v.ok.map(b => if (b) "1" else ".").mkString))
				
		println("RATIO VALIDATIONS: ")
		println("          \t"+("".padTo(c.g.nistart, '-').padTo(c.g.niend, '#').padTo(c.g.iend - c.g.istart + 1, '-')))
		Ratios.iterate(c.fragmentValids.rValids.length, (ri, ui, di) => {
			println(aIn.cg.chromatograms(ui).q3.toInt + " / " + 
					aIn.cg.chromatograms(di).q3.toInt + " \t" +
					c.fragmentValids.rValids(ui)(di).ok.map(b => if (b) "1" else ".").mkString +
					" \t%.1e".format(c.g.pvalueTable(ui)(di))
				)
		})
		
		fragGraph.repaint
		
		
		//if (g.estimates == null) g.estimateAndIntegrate(state, aIn.cg.chromatograms.map(_.intensities).toArray)
		showEstimates(c)
	}
	
	
	
	def showEstimates(c:Carrier) = {
		val i0 = c.g.istart - 8
		val in = c.g.iend + 8
		
		val origTimes = aIn.cg.chromatograms.head.times.slice(i0, in)
		val estTimes = origTimes.slice(8+c.g.nistart, 8+c.g.niend)
		
		var origCurves = aIn.cg.chromatograms.map(c => 
			new Curve2(
				origTimes,
				c.intensities.slice(i0, in),
				origTimes.map(_ => false),
				null,
				Color.LIGHT_GRAY
			))
		var estCurves = aIn.cg.chromatograms.zip(c.fragmentEstimation.estimates).map(tu => 
			new Curve2(
				estTimes,
				tu._2,
				estTimes.map(_ => false),
				null,
				Color.RED
			))
			
		estimateGraph.title = c.g+" estimates"
		estimateGraph.setCurves(origCurves ++ estCurves)
		estimateGraph.repaint
	}
	
	
	
	def showRef(aIn:AnubisInput) = {
		val xs = Array(-2.0, -1.0, 0.0, 1.0, 2.0)
		var cs = aIn.refRatios.filter(_.transitionId1 == 0).map(r => {
			new Curve2(
				xs, 
				xs.map(x => (4 - x*x)/r.mean), 
				xs.map(d => false), 
				"q3:"+aIn.cg.chromatograms(r.transitionId2).q3
			)
		})
		refGraph.setCurves(new Curve2(
				xs, 
				xs.map(x => (4 - x*x)), 
				xs.map(d => false), 
				"q3:"+aIn.cg.chromatograms(0).q3
			) +: cs)
		
		refGraph.setZoom(-2.5, 2.5)
		refGraph.repaint
	}
	

	
	
	def printRatioStats(cg:XChromatogramGroup):(Double, Array[Array[Double]]) = {
		
		val cL 		= cg.chromatograms.length
		var rL		= pepQuant.ratios.length
		var pearson = new PearsonsCorrelation
		var ratios	= pepQuant.ratios
		
		
		
		var max 	= Matrix.get2d[Double](cL)
		var median 	= Matrix.get2d[Double](cL)
		var v 		= Matrix.get2d[Double](cL)
		var vSum 	= 0.0
		
		Ratios.iterate(cL, (ri, ui, di) => {
			var r = ratios.getRatio(ri)
						.filter(d => d > 0.0 && d < Double.PositiveInfinity)
			max(ui)(di) = r.max
			median(ui)(di) = StatUtils.percentile(r, 50)
			v(ui)(di) = StatUtils.variance(r)
			vSum += v(ui)(di)
		})
		vSum /= 45
		
		/*
		println("RATIO MAX ("+cL+"):")
		printMatrix(max)
		println("RATIO MEDIAN ("+cL+"):")
		printMatrix(median)
		println("RATIO VARIANCES ("+cL+"):")
		printMatrix(v)
		*/
			
			
		
		var corr 	= Matrix.get2d[Double](rL)
		var corr3 	= Matrix.get2d[Double](rL)
		
		Ratios.iterate(rL, (ri, ui, di) => {
			var rr = ratios.getRatio(ui).zip(ratios.getRatio(di))
						.filter(t => 
									t._1 > 0.0 && t._1 < Double.PositiveInfinity
								&&	t._2 > 0.0 && t._2 < Double.PositiveInfinity
							).unzip
			corr(ui)(di) = pearson.correlation(rr._1.toArray, rr._2.toArray)
			corr3(ui)(di) = pearson.correlation(rr._1.map(math.log _).toArray, rr._2.map(math.log _).toArray)
		})
		
		/*
		println("RATIO CORRELATIONS ("+rL+"):")
		printMatrix(corr)
		println("RATIO CORRELATIONS LOGGED ("+rL+"):")
		printMatrix(corr3)
		*/
		
		return (vSum, corr)
	}
	
	
	
	def printMatrix(m:Array[Array[Double]]):Unit = {
		var s = m.map(_.map(d => if (d <= 0.0) "." else math.log10(math.abs(d)).toInt.toString))
		var width = s.map(_.map(_.length).max).max
		for (i <- 0 until s.length)
			println(s(i).map(x => "".padTo(width - x.length, ' ') + x).mkString(" "))
	}
	
	
	
	
	def clearDenoise = {
		denoiseGraph.clear
		denoiseGraph.repaint
	}
	
	
	
	
	def status(str:String) = statusL.text = str
}