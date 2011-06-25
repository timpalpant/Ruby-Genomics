require 'fixed_precision'
require 'genomic_index_error'
require 'unix_file_utils'
require 'enumerator'

##
# Represents a contiguous block of genomic data values
# Returned when querying a WigFile / SpotArray
#
##
class Contig
  include Enumerable
  attr_accessor :chr

  def initialize(chr = 'unknown')
    @chr = chr
    @data = Hash.new
  end
  
  def each
    (start..stop).each { |bp| yield get(bp) }
  end
  
  ##
  # ACCESS METHODS
  ##
  
  def set(base, value)
    @data[base] = value
  end
  
  def get(base)
    @data[base]
  end
  
  # Analogous to Array#[]
  # Can take a single integer (base pair), or a range, or two integers (slice)
  def [](*args)
    if args.length == 1
      arg = args.first
      if arg.is_a?(Integer)
        raise ContigError, "Contig does not contain data for the base (#{arg})" unless include?(arg)
        return get(arg)
      elsif arg.is_a?(Range)
        raise ContigError, "Contig does not contain data for the range #{arg}" unless include?(arg.min, arg.max)
        return arg.map { |base| get(base) }
      else
        raise ContigError, "Invalid type of argument passed to Contig (#{arg.class})"
      end
    elsif args.length == 2
      start = args[0]
      length = args[1]
      stop = start + length - 1
      raise ContigError, "Invalid type of arguments passed to Contig (#{start.class}, #{length.class})" unless start.is_a?(Integer) and length.is_a?(Integer)
      raise ContigError, "Contig does not contain data for the range (#{start}..#{start+length-1})" unless include?(start, start+length-1)
      
      return (start..stop).map { |base| get(base) }
    else
      raise ContigError, "Invalid number of arguments passed to Contig (1 or 2 args accepted, #{args.length} passed!)"
    end
  end
  
  ##
  # PROPERTY METHODS
  ##
  
  # The first base pair with data
  def start
    @data.keys.min
  end
  
  # The last base pair with data
  def stop
    @data.keys.max
  end
  
  # The number of base pairs of data
  # NOTE: some bases could be missing data if the block is not truly contiguous
  def length
    stop - start + 1
  end
  
  # If this Contig contains data for the specified base pair / range
  def include?(low, high = self.length-1)
    low >= start and high <= stop
  end
  
  # Alias for #include
  def cover?(low, high = self.length-1)
    include?(low, high)
  end
  
  # Get a subsequence of data as an Array
  # Alias for #[], but also allows Crick querying
  def bases(from, to)
    low = [from, to].min
    high = [from, to].max
    
    # Get the values
    values = self[low..high]
    
    # Allow crick querying
    values.reverse! if from > to
    return values
  end
  
  # Alias for query
  def bases(from, to)
    query(from, to)
  end
  
  # Alias for query
  def subseq(low, high)
    query(low, high)
  end

  ##
  # OUTPUT METHODS
  ##
  
  # Convert this Contig into an Array of values
  def to_a
    (start..stop).map { |bp| get(bp) }
  end
  
  # Output this Contig as a variableStep Wiggle block
  def to_variable_step
    str = StringIO.new
    str << "variableStep chrom=#{@chr} span=1"
    
    @data.each do |bp,value|
      str << "\n#{bp}\t#{value}"
    end
    
    return str.string
  end
  
  # Output this Contig as a fixedStep Wiggle block
  def to_fixed_step
    str = StringIO.new
    str << "fixedStep chrom=#{@chr}"
    str << " start=#{start} step=1 span=1"

    self.each do |value|
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
  def to_contig(chr = 'unknown', start = 1)
    contig = Contig.new(chr)
    self.each_with_index { |value,i| contig.set(start+i, value) }
    return contig
  end
end


# Raised if something goes wrong with a Contig
class ContigError < StandardError
end
