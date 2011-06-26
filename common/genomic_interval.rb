##
# "Abstract" class for all types of genomic coordinates
# Subclasses: Spot, Read
##
class GenomicInterval
  attr_accessor :chr, :start, :stop
  
  def initialize(chr = nil, start = nil, stop = nil)
    @chr = chr
    @start = start
    @stop = stop
  end
  
  # Output this genomic interval as a string (the format is suitable
  # for UCSC or samtools)
  def to_s
    "#{@chr}:#{start}-#{stop}"
  end
  
  def to_bed
    "#{@chr}\t#{low}\t#{high}\t.\t.\t#{strand}"
  end
  
  def to_bedgraph
    "#{@chr}\t#{low}\t#{high}"
  end
  
  def to_gff
    "#{@chr}\tSpotArray\tfeature\t#{low}\t#{high}\t0\t.\t.\tprobe_id=none;count=1"
  end
  
  def center
    (@start + @stop) / 2 if @start and @stop
  end
  
  def length
    (@stop - @start).abs + 1 if @start and @stop
  end
  
  # Whether this spot includes (encompasses) a given locus
  def include?(base)
    low <= base and high >= base
  end
 
  # The minimum chromosomal coordinate, regardless of strand
  def low
    # Cache for performance
    @low = [@start, @stop].min if @low.nil?
    
    return @low
  end
 
  # The maximum chromosomal coordinate, regardless of strand
  def high
    # Cache for performance
    @high = [@start, @stop].max if @high.nil?
    
    return @high
  end
 
  # If this spot is oriented on the plus strand (Watson)
  def watson?
    @stop >= @start
  end
 
  def crick?
    not watson?
  end
  
  def strand
    watson? ? '+' : '-'
  end
  
  # TODO: Other conditions for being valid?
  def valid?
    @start > 0 and @stop > 0
  end
end