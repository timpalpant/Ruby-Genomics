require 'genomic_interval'

##
# A genomic interval with an associated value (such as a BedGraph entry)
# i.e. for microarrays
##
class Spot < GenomicInterval
  attr_accessor :id, :value
  
  def initialize(id = nil)
    @id = id
  end
  
  def to_s
    "Spot: #{@id},#{chr},#{@start},#{@stop},#{@value}"
  end
  
  def to_bed
    "#{chr}\t#{@start}\t#{@stop}\t#{@id}\t#{@value}"
  end
  
  def to_bedGraph
    "#{chr}\t#{@start}\t#{@stop}\t#{@value}"
  end
  
  def to_gff
    "#{chr}\tSpotArray\tfeature\t#{@start}\t#{@stop}\t#{@value}\t.\t.\tprobe_id=#{@id};count=1"
  end
end
