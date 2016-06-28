package se.lth.immun.diana

import akka.actor.Actor

import se.lth.immun.graphs.LineGraph
import se.lth.immun.graphs.util._

import java.io.File
import java.io.IOException
import java.awt.Color
import java.awt.Dimension
import java.awt.image.BufferedImage
import javax.imageio.ImageIO

object AssayReportActor {
	case class WriteReport(da:DianaAssay, res:Seq[Carrier], aIn:DianaInput, iDIn:Option[DianaInput])
}

class AssayReportActor(params:DianaScorerParams) extends Actor {

	import AssayReportActor._
	
	
	def receive = {
		case WriteReport(da, res, aIn, iDIn) =>
			
			val c = res.minBy(_.fragmentRankPcsRatioProb)
			val i0 = c.g.istart - 16
			val in = c.g.iend + 16
			
			val fg = fragGraph(da, c, aIn, i0, in)
			val ig = iDIn.map(x => isoGraph(da, c, x, i0, in))
			
			val image = new BufferedImage(400, 400, BufferedImage.TYPE_INT_RGB)
			val g = image.createGraphics
			g.setColor(Color.WHITE)
			g.fillRect(50, 50, 50, 80)
			
			fg.renderer.setup(fg.xAxis, fg.yAxis, fg.style, new Size(400, 200))
			fg.render(g, fg.renderer)
		
			ig match {
				case Some(lg) =>
					g.translate(0, 200)
					lg.renderer.setup(lg.xAxis, lg.yAxis, lg.style, new Size(400, 200))
					lg.render(g, lg.renderer)
				case None => {}
			}
			
			g.dispose
			try { 
			    ImageIO.write(image, "png", new File("qc/"+da.pepCompRef.replace("/","z").replace(':','_')+".png")) 
			} catch {
				case ioe:IOException =>
			    	ioe.printStackTrace
			}
			
	}
	
	
	def fragGraph(da:DianaAssay, c:Carrier, aIn:DianaInput, i0:Int, in:Int) = {
		val fg 	= new LineGraph
		fg.preferredSize 			= new Dimension(200, 200)
		fg.style.annotColBackground = new Color(0.5f, 0.5f, 0.5f, 0.5f)
		
		val origTimes = aIn.cg.chromatograms.head.times.slice(i0, in)
		//val estTimes = origTimes.slice(8+c.g.nistart, 8+c.g.niend)
		
		var origCurves = aIn.cg.chromatograms.map(c => 
			new Curve2(
				origTimes,
				c.intensities.slice(i0, in),
				origTimes.map(_ => false),
				c.id.split("/").last
			))
		/*var estCurves = aIn.cg.chromatograms.zip(c.fragmentEstimation.estimates).map(tu => 
			new Curve2(
				estTimes,
				tu._2,
				estTimes.map(_ => false),
				tu._1.id.split("/").last
			))
			*/
		fg.title = "fragments"
		fg.setCurves(origCurves) // ++ estCurves)
		fg.addAnnotation(
				new HeightBoxAnnotation(
					origTimes(math.max(i0, math.min(c.g.istart, in))-i0), 
					c.fragmentEstimation.estimates.map(_.max).max, 
					origTimes(math.max(i0, math.min(c.g.iend, in))-i0), 
					Annotation.BACKGROUND
				))
		fg.repaint
		fg
	}
	
	
	def isoGraph(da:DianaAssay, c:Carrier, iDIn:DianaInput, i0:Int, in:Int) = {
		val fg 	= new LineGraph
		fg.preferredSize 			= new Dimension(200, 200)
		fg.style.annotColBackground = new Color(0.5f, 0.5f, 0.5f, 0.5f)
		
		val origTimes = iDIn.cg.chromatograms.head.times.slice(i0, in)
		val estTimes = origTimes.slice(8+c.g.nistart, 8+c.g.niend)
		
		var origCurves = iDIn.cg.chromatograms.map(c => 
			new Curve2(
				origTimes,
				c.intensities.slice(i0, in),
				origTimes.map(_ => false),
				c.id.split(" ", 2).last
			))
			
		fg.title = "prec isotopes"
		fg.setCurves(origCurves)
		fg.repaint
		fg
	}
}
