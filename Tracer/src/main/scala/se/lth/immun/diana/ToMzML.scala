package se.lth.immun.diana

import java.util.Arrays

import se.lth.immun.base64.Base64
import se.lth.immun.mzml.BinaryDataArray
import se.lth.immun.mzml.CvParam
import se.lth.immun.mzml.ghost.GhostBinaryDataArray
import se.lth.immun.mzml.ghost.Ghost
import se.lth.immun.mzml.ghost.GhostException
import ms.numpress.MSNumpress

object ToMzML {

	def toBinaryDataArray(
			bytes:Iterable[Byte],
			dd:GhostBinaryDataArray.DataDef
	):BinaryDataArray = {
		
		var b = new BinaryDataArray
		annotate(b, dd)
		b.binary 		= Base64.decoder.encodeToString(bytes.toArray)
		b.encodedLength = b.binary.length
		
		b
	}
	
	def annotate(b:BinaryDataArray, dd:GhostBinaryDataArray.DataDef) = {
		import Ghost._
		import GhostBinaryDataArray._
		import MSNumpress._
		
		dd.dataType match {
			case Time() =>
				b.cvParams += new CvParam {
					cvRef = "MS"
					accession = TIME_ARRAY_ACC
					name = "time array"
					unitCvRef = Some("UO")
					unitAccession = Some("UO:0000010")
					unitName = Some("second")
				}
			case Intensity() =>
				b.cvParams += new CvParam {
					cvRef = "MS"
					accession = INTENSITY_ARRAY_ACC
					name = "intensity array"
					unitCvRef = Some("MS")
					unitAccession = Some("MS:1000131")
					unitName = Some("number of counts")
				}
			case MZ() =>
				b.cvParams += new CvParam {
					cvRef = "MS"
					accession = INTENSITY_ARRAY_ACC
					name = "intensity array"
					unitCvRef = Some("MS")
					unitAccession = Some("MS:1000131")
					unitName = Some("number of counts")
				}
			case x =>
				throw new GhostException("Cannot write binary data of unknown type "+x)
		}
		
		if (dd.zlibCompression)
			b.cvParams += new CvParam {
				cvRef = "MS"
				accession = ZLIB_COMPRESSION_ACC
				name = "zlib compression"
			}
		
		if (dd.doublePrecision)
			b.cvParams += new CvParam {
				cvRef = "MS"
				accession = BIT_64_ACC
				name = "64-bit float"
			}
		else
			b.cvParams += new CvParam {
				cvRef = "MS"
				accession = BIT_32_ACC
				name = "32-bit float"
			}
		
		val numpressed = dd.numCompression match {
			case ACC_NUMPRESS_LINEAR => 
				b.cvParams += new CvParam {
					cvRef = "MS"
					accession = ACC_NUMPRESS_LINEAR
					name = "MS-Numpress linear prediction compression"
				}
				true
								
			case ACC_NUMPRESS_PIC => 
				b.cvParams += new CvParam {
					cvRef = "MS"
					accession = ACC_NUMPRESS_PIC
					name = "MS-Numpress positive integer compression"
				}
				true
				
			case ACC_NUMPRESS_SLOF => 
				b.cvParams += new CvParam {
					cvRef = "MS"
					accession = ACC_NUMPRESS_SLOF
					name = "MS-Numpress short logged float compression"
				}
				true
				
			case _ => false
		}
		
		if (!numpressed && !dd.zlibCompression)
			b.cvParams += new CvParam {
				cvRef = "MS"
				accession = NO_COMPRESSION_ACC
				name = "no compression"
			}
	}
	
}