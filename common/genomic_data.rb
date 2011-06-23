##
# Abstract class for all arrays of genomic information
# stored by chromosome
# Subclasses: Genome, Assembly
#
# The structure is a Hash of chromosomes (keys), with an array of 
# spots/reads in each value
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
