require 'fixed_precision'
require 'sparse_array'

##
# Represents a contiguous block of genomic data values
# Returned when querying a WigFile / SpotArray
#
##
class Contig < SparseArray
  attr_accessor :chr

  def initialize(chr = 'unknown')
    super()
    @chr = chr
  end

  ##
  # OUTPUT METHODS
  ##
  
  # Output this Contig as a variableStep Wiggle block
  def to_variable_step(res = 1)
    str = StringIO.new
    str << "variableStep chrom=#{@chr} span=#{res}"
    
    (start..stop).step(res) do |bp|
      value = get(bp)
      str << "\n#{bp}\t#{value.to_s(5)}" if value
    end
    
    return str.string
  end
  
  # Output this Contig as a fixedStep Wiggle block
  def to_fixed_step(res = 1)
    str = StringIO.new
    str << "fixedStep chrom=#{@chr}"
    str << " start=#{start} step=#{res} span=#{res}"

    (start..stop).step(res) do |bp|
      value = get(bp)
      if value 
        str << "\n" << value.to_s(5)
      else
        str << "\nNaN"
      end
    end

    return str.string
  end
end


# For converting an Array to a Contig
class Array
  def to_contig(chr = 'unknown', start = 1, step = 1, span = 1)
    contig = Contig.new(chr)
    
    bp = start
    self.each do |value|
      (bp..bp+span-1).each { |base| contig.set(base, value) }
      bp += step
    end
    
    return contig
  end
end


# Raised if something goes wrong with a Contig
class ContigError < StandardError
end
