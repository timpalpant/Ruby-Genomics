#
#  sparse_array.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

#  Stores info like an Array, but internally may manage it differently to be more memory-efficient
#  Provides a layer of abstraction between the storage (Array, Hash, etc.) and the downstream users (Contig)
#  that allows us to flexibly change implementation details to optimize performance without affecting other classes

require 'enumerator'
require 'stats'

class SparseArray
  include Enumerable
  
  def initialize
    @data = Hash.new
  end
  
  def each
    (start..stop).each { |i| yield get(i) }
  end
  
  ##
  # ACCESS METHODS
  ##
  
  # Store a value for a given index
  def set(i, value)
    # Don't store nil values
    return if value.nil?
    
    raise SparseArrayError, "Invalid key (#{i.class})! Indexes must be Integers" unless i.is_a?(Integer)
    raise SparseArrayError, "Invalid key (#{i})! Cannot set value of index < 0" if i < 0
    
    # Store the value
    @data[i] = value
    
    # Update the cached min/max if it has changed
    @min = i if @min.nil? or i < @min
    @max = i if @max.nil? or i > @max
  end
  
  # Alias for #set
  def []=(i,value)
    set(i, value)
  end
  
  # Get an object for an index, or return nil if we don't have data
  def get(i)
    raise SparseArrayError, "SparseArrays do not support negative indexing" if i < 0
    @data[i]
  end
  
  # Alias for #get
  def fetch(i)
    get(i)
  end
  
  # Analogous to Array#[]
  # Can take a single integer (base pair), or a range, or two integers (slice)
  def [](*args)
    if args.length == 1
      arg = args.first
      if arg.is_a?(Integer)
        raise SparseArrayError, "SparseArrays do not support negative indexing" if arg < 0
        raise SparseArrayError, "SparseArray does not contain data for the base (#{arg})" unless cover?(arg)
        return get(arg)
      elsif arg.is_a?(Range)
        raise SparseArrayError, "SparseArrays do not support negative indexing" if arg.min < 0 or arg.max < 0
        raise SparseArrayError, "SparseArray does not contain data for the range #{arg}" unless cover?(arg.min, arg.max)
        return arg.map { |base| get(base) }
      else
        raise SparseArrayError, "Invalid type of argument passed to SparseArray (#{arg.class})"
      end
    elsif args.length == 2
      start = args[0]
      length = args[1]
      raise SparseArrayError, "Invalid type of arguments passed to SparseArray (#{start.class}, #{length.class})" unless start.is_a?(Integer) and length.is_a?(Integer)
      raise SparseArrayError, "SparseArrays do not support negative indexing" if start < 0

      return Array.new if length == 0
      
      stop = start + length - 1
      raise SparseArrayError, "SparseArray does not contain data for the range (#{start}..#{start+length-1})" unless cover?(start, start+length-1)
      
      return (start..stop).map { |base| get(base) }
    else
      raise SparseArrayError, "Invalid number of arguments passed to SparseArray (1 or 2 args accepted, #{args.length} passed!)"
    end
  end
  
  # Get all of the indices as an Array
  def indices
    @data.keys
  end
  
  # Alias for indices
  def keys
    indices
  end
  
  # Get all of the values as an Array
  def values
    (start..stop).map { |i| get(i) }
  end
  
  ##
  # PROPERTY METHODS
  ##
  
  # The lowest index with data
  def start
    @min
  end
  
  # Alias for start
  def low
    start
  end
  
  # The highest index with data
  def stop
    @max
  end
  
  # Alias for stop
  def high
    stop
  end
  
  # The number of indices from start..stop
  def length
    stop - start + 1 unless coverage == 0
  end
  
  # Number of indices with data
  def coverage
    @data.length
  end
  
  # A measure of how sparse this array is
  # lower values = more sparse
  def sparsity
    coverage.to_f / length
  end
  
  # If this SparseArray contains data encompassing the specified index / range
  def include?(low, high = nil)
    return @data.include?(low) if high.nil?
    (low..high).each { |i| return false if get(i).nil? }
    return true
  end
  
  # If the specified base / range falls with the start-stop of this SparseArray
  def cover?(low, high = @max)
    return false if coverage == 0 or low.nil?
    @min <= low and low <= @max and high >= @min and high <= @max
  end
  
  # Get a subsequence of data as an Array
  # Alias for #[], but also allows reverse querying
  def query(from, to)
    low = [from, to].min
    high = [from, to].max
    
    # Get the values
    values = self[low..high]
    
    # Allow crick querying
    values.reverse! if from > to
    return values
  end
  
  ##
  # MATH
  ##
  
  def +(other)
    binary_op(other) { |value1,value2| value1+value2 }
  end

  def -(other)
    binary_op(other) { |value1,value2| value1-value2 }
  end
  
  def *(other)
    binary_op(other) { |value1,value2| value1*value2 }
  end

  def /(other)
    binary_op(other) { |value1,value2| value1/value2 }
  end

  def sum
    @data.values.sum
  end

  # Alias for #sum
  def total
    sum
  end

  def mean
    @data.values.mean
  end

  def stdev
    @data.values.stdev
  end

  def median
    @data.values.median
  end

  ##
  # OUTPUT METHODS
  ##
  
  # Convert this SparseArray into an Array of values
  # Alias for #values
  def to_a
    values
  end
  
  # Convert this SparseArray into a Hash of key-value pairs
  def to_hash
    @data
  end

  ##
  # HELPER METHODS
  ##

  private

  def binary_op(other)
    result = SparseArray.new
    
    if other.is_a?(Integer) or other.is_a?(Float)
      @data.each { |i,value| result.set(i, yield(value, other)) }
    elsif other.is_a?(SparseArray)
      min_i = [self.start, other.start].min
      max_i = [self.stop, other.stop].max
      
      (min_i..max_i).each do |i|
        result.set(i, yield(self[i], other[i])) if self.include?(i) and other.include?(i)
      end
    else
      raise SparseArrayError, "Invalid argument of type #{other.class}!"
    end
    
    return result
  end
end

class SparseArrayError < StandardError
end