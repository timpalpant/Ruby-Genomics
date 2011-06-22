require 'bio'
require 'genomic_interval'

##
# A genomic interval with an associated sequence (read)
# i.e. for high-throughput, short-read sequencing
##
class Read < GenomicInterval
  attr_accessor :seq, :qual
	
	def to_s
		"Read: #{chr},#{@start},#{@stop}"
	end
	
	def to_bed
		"#{chr}\t#{@start}\t#{@stop}"
	end
  
  def to_bedgraph
		to_bed
	end
  
  def to_sam
  	raise "Not yet implemented"
  end
end
