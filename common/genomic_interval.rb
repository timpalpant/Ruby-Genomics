##
# "Abstract" class for all types of genomic coordinates
# Subclasses: Spot, Nucleosome, Read, Locus, etc.
##
class GenomicInterval
	attr_accessor :start, :stop
	
	def initialize(start = nil, stop = nil)
		@start = start
		@stop = stop
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
	
	# TODO: Other conditions for being valid?
	def valid?
		@start > 0 and @stop > 0
	end
end