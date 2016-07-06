#!/usr/bin/env python2
"""
Copy from Allan Haldane.
"""
from __future__ import with_statement
from scipy import *
import numpy
import sys, types
import pylab
import matplotlib.cm
from matplotlib.patches import PathPatch
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
from matplotlib.font_manager import FontProperties

#helper class that opens a file if it is not already opened
class Opener:
	def __init__(self, fileobj, rw="rt"):
		self.fileobj = fileobj
		self.rw = rw
		self.f = None

	def __enter__(self):
		if isinstance(self.fileobj, types.StringType):
			self.f = open(self.fileobj, self.rw)
			return self.f
		elif hasattr(self.fileobj, 'read') or hasattr(self.fileobj, 'write'):
			if self.rw != self.fileobj.mode: #error? XXX
				raise Exception("File is already open, but in wrong mode")
			return self.fileobj

	def __exit__(self, e, fileobj, t):
		if self.f != None:
			self.f.close()
		return False

def loadSites(fn, alphabet): #optimized for fast loading, assumes ASCII
	with Opener(fn) as f:
		seqs = [l.strip() for l in f if not l.startswith('#') and l != '']

	nSeqs, seqLen = len(seqs), len(seqs[0])
	if any([len(s) != seqLen for s in seqs]): #check that file is OK
		raise Exception("Error: Sequences have different lengths")

	nucNums = -ones(256, int) #nucNums is a map from ascii to base number
	nucNums[frombuffer(alphabet, uint8)] = arange(len(alphabet))

	charnums = frombuffer("".join(seqs), uint8).reshape((nSeqs, seqLen))
	seqTable = nucNums[charnums]

	badBases = seqTable < 0
	if any(badBases):
		bad = " ".join(chr(i) for i in set(charnums[badBases]))
		raise Exception("Unknown characters: {0}".format(bad))

	return seqTable

def getCounts(seqs, nBases):
	nSeq, seqLen = seqs.shape
	bins = arange(nBases+1, dtype='int')
	counts = zeros((seqLen, nBases), dtype='int')
	for i in range(seqLen):
		counts[i,:] = histogram(seqs[:,i], bins)[0]
	return counts # index as [pos, res]

def getFreqs(seq, nBases):
	return getCounts(seq, nBases).astype('float')/seq.shape[0]

def showFreqs(f, alphabet="ACGT", nSeqs=inf, fig=None):
	#according to logo standard
	seqLen, nBases = f.shape
	if nBases != len(alphabet):
		raise Exception("Need %d alphabet" % len(alphabet))

	if nBases == 4:
		colors = ["#00cc00", "#0000cc", "#ffb300", "#cc0000"]
	else:
		colors = matplotlib.cm.jet(linspace(0,1,nBases))

	fontprop = FontProperties(family="monospace", stretch="condensed",
							  weight="bold", size="medium")

	#disable warnings for a few lines, because of the log2
	errsettings = numpy.seterr(all='ignore')

	#see the scheider papers on logos for description
	#here we use the 'approximate' method, good for nseq > 50.
	ecorr = (nBases-1)/(2*log(2)*nSeqs)
	f = f/sum(f,1)[:,newaxis]
	R = log2(nBases) - (-sum(nan_to_num(f*log2(f)),1)+ecorr)
	heights = f#*R[:,newaxis]

	numpy.seterr(**errsettings)

	if fig == None:
		fig = pylab.figure(figsize=(8,3))

	def addLetter(char, x, y, height, color):
		text = TextPath((0,0), char, size=1, prop=fontprop)
		bbox = text.get_extents()
		(tx0,ty0), tw, th = bbox.min, bbox.width, bbox.height
		trans = Affine2D.identity() \
						.translate(-tx0, -ty0) \
						.scale(0.9/tw, height/th) \
						.translate(x+0.5, y)
		t = trans.transform_path(text)
		pylab.gca().add_patch(PathPatch(t, fc=color, ec='none'))

	for n,h in enumerate(heights):
		letters = [(alphabet[b], h[b], colors[b]) for b in argsort(h)]

		ypos = 0
		for char, height, color in letters:
			addLetter(char, n, ypos, height, color)
			ypos += height

	pylab.xticks(arange(seqLen)+1)
	pylab.ylabel('Bits')
	pylab.xlim(0.5,seqLen+0.5)
	#pylab.ylim(0,max(sum(heights,1)))
	#pylab.ylim(0,log2(nBases))
        pylab.ylim(0, 1)
	return fig

def showLogo():
	#try to load either freq table or em:
	try:
		if len(sys.argv) == 4:
			e = float(sys.argv[3])
			f = exp(-e*loadtxt(sys.argv[1]))
		else:
			f = loadtxt(sys.argv[1])
		alphabet = sys.argv[2]
		showFreqs(f, alphabet)
		pylab.show()
		return
	except Exception as e:
		pass

	#if that failed, try loading sequences:
	if len(sys.argv) > 2:
		trialalphabet = [sys.argv[2]]
	else:
		trialalphabet = []
	trialalphabet.extend(["ACGT", "ACDEFGHIKLMNPQRSTVWY"])

	if len(sys.argv) > 1:
		for alphabet in trialalphabet:
			try:
				seq = loadSites(sys.argv[1], alphabet)
			except Exception as e:
				continue
			f = getFreqs(seq, len(alphabet))
			showFreqs(f, alphabet, seq.shape[0])
			pylab.show()
			return
		print 'Error: Could not detect alphabet in {}'.format(sys.argv[1])
		return

	#otherwise print usage
	print "Usage: one of:"
	print "logos.py sequence_file     (autodetects DNA or protein sequences)"
	print "logos.py sequence_file alphabet"
	print "logos.py freq_file alphabet"
	print "logos.py em_file alphabet B"

def showEnergies(em, alphabet="ACGT", nSeqs=inf, fig=None):
	#according to logo standard
	seqLen, nBases = em.shape
	if nBases != len(alphabet):
		raise Exception("Need %d alphabet" % len(alphabet))

	if nBases == 4:
		colors = ["#00cc00", "#0000cc", "#ffb300", "#cc0000"]
	else:
		colors = matplotlib.cm.jet(arange(1,nBases+1, dtype=float)/(nBases+1))

	fontprop = FontProperties(family="monospace", stretch="condensed",
							  weight="bold", size="medium")
	em = em-mean(em,1)[:,newaxis]
	print em
	sem = argsort(em,0)

	if fig == None:
		fig = pylab.figure(figsize=(8,3))

	def addletter(let, x, y, height, color, alpha=1):
		text = TextPath((0,0),let,size=1, prop=fontprop)
		bbox = text.get_extents()
		tx0,ty0 = bbox.min
		tw,th = bbox.width, bbox.height
		trans = Affine2D.identity() \
						.translate(-tx0,-ty0) \
						.scale(0.9/tw,height/th) \
						.translate(x,y)
		t = trans.transform_path(text)
		pylab.gca().add_patch(PathPatch(t, fc=color, ec='none',alpha=alpha))

	for n,row in enumerate(em):
		pos = arange(nBases)[row>0][argsort(row[row>0])]
		neg = arange(nBases)[row<0][argsort(row[row<0])[::-1]]

		hpos = 0
		for i in pos:
			addletter(alphabet[i], n+0.5, hpos, row[i], colors[i])
			hpos += row[i]
		hpos = 0
		for i in neg:
			addletter(alphabet[i], n+0.5, hpos+row[i], -row[i], colors[i],alpha=0.5)
			hpos += row[i]

	pylab.hlines(0, 0, seqLen+0.5, colors='k')
	pylab.xticks(arange(seqLen)+1)
	#pylab.ylabel('Bits')
	pylab.xlim(0.5,seqLen+0.5)
	#pylab.ylim(0,max(sum(heights,1)))
	#pylab.ylim(0,log2(nBases))
	return fig

def showAffinity():
	f = loadtxt(sys.argv[1])
	if len(sys.argv) == 2:
		if f.shape[1] == 4:
			alphabet = 'ACGT'
		else:
			raise Exception("Unknown alphabet!")
	else:
		alphabet = sys.argv[2]
	fig = showEnergies(f, alphabet)
	fig.subplots_adjust(left=0.1,right=0.9)
	pylab.show()

def main():
	showLogo()
	#pylab.savefig('filename')
	#showAffinity()

if __name__ == '__main__':
	main()

