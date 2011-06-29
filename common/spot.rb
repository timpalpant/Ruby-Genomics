require 'genomic_interval'
require 'stringio'

##
# A genomic interval with an associated value (such as a BedGraph entry)
# i.e. for microarrays
##
class Spot < GenomicInterval
  attr_accessor :id, :value
  
  def initialize(chr = nil, start = nil, stop = nil, id = nil, value = nil)
    super(chr, start, stop)
    @id = id
    @value = value
  end
  
  def to_s
    "Spot: #{@id},#{chr},#{@start},#{@stop},#{@value}"
  end
  
  def to_bed
    value = @value ? @value : '.'
    id = @id ? @id : '.'
    "#{@chr}\t#{low-1}\t#{high}\t#{id}\t#{value}\t#{strand}"
  end
  
  def to_bedgraph
    bedgraph = StringIO.new
    bedgraph << "#{chr}\t#{low-1}\t#{high}"
    bedgraph << "\t#{@value}" if @value    
    bedgraph.string
  end
  
  def to_gff
    id = @id ? @id : 'no_id'
    value = @value ? @value : '.'
    "#{@chr}\tSpotArray\tfeature\t#{low}\t#{high}\t#{value}\t#{strand}\t.\tprobe_id=#{id};count=1"
  end
end
