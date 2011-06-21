##
# Abstract class for all arrays of genomic information
# stored by chromosome
# Subclasses: SpotArray, NukeCalls, SAM, Bed, BedGraph, etc.
#
# The structure is a Hash of chromosomes, with an array of 
# spots/reads/entries under each Hash
##
class GenomicData < Hash
	# Return all of the chromosomes in this genome
	def chromosomes
		self.keys
	end
	
	# Return a specific chromosome
	def chr(chr_id)
		self[chr_id]
	end
end
