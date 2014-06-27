package se.lth.immun.diana

import se.lth.immun.unimod._
import se.lth.immun.chem._

object PeptideParser {

	trait PepParseResult
	case class UniModPeptide(p:Peptide) extends PepParseResult
	case class Unparsable(seq:String) extends PepParseResult
	
	def parseSequence(seq:String):PepParseResult = {
		try {
			return UniModPeptide(UniMod.parseUniModSequence(seq))
		} catch {
			case e:Throwable =>
				return Unparsable(seq)
		}	
	}
}