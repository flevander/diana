package se.lth.immun.diana

import scala.collection.mutable.ArrayBuilder
import scala.collection.mutable.WrappedArray

object TraceBuilder {
	def apply() = new TraceBuilder()
}

class TraceBuilder(val sizeIncrement: Double = 2.0) extends ArrayBuilder[Double] {

	protected var elems = new Array[Double](0)
	protected var capacity: Int = 0
	protected var size: Int = 0

	private def mkArray(size: Int): Array[Double] = {
		val newelems = new Array[Double](size)
		if (this.size > 0) Array.copy(elems, 0, newelems, 0, this.size)
		newelems
	}

	private def resize(size: Int) {
		elems = mkArray(size)
		capacity = size
	}

	override def sizeHint(size: Int) {
		if (capacity < size) resize(size)
	}

	protected def ensureSize(size: Int) {
		if (capacity < size || capacity == 0) {
			var newsize = if (capacity == 0) 16 else capacity * sizeIncrement
			while (newsize < size) newsize *= sizeIncrement
			resize(newsize.toInt + 1)
		}
	}

	def +=(elem: Double): this.type = {
		ensureSize(size + 1)
		elems(size) = elem
		size += 1
		this
	}

	override def ++=(xs: TraversableOnce[Double]): this.type = xs match {
		case xs: WrappedArray.ofDouble =>
			ensureSize(this.size + xs.length)
			Array.copy(xs.array, 0, elems, this.size, xs.length)
			size += xs.length
			this
		case _ =>
			super.++=(xs)
	}

	def clear() {
		size = 0
	}

	def result() = {
		if (capacity != 0 && capacity == size) elems
		else mkArray(size)
	}

	def nocopyResult = elems.view(0, size)

	override def equals(other: Any): Boolean = other match {
		case x: TraceBuilder => (size == x.size) && (elems == x.elems)
		case _ => false
	}

	override def toString = "TraceBuilder"
}