require 'assembly'
require 'ucsc_tools'
require 'stats'
require 'enumerator'

##
# Get genomic sequences from 2bit files
# Loads the sequences as Bio::Sequence::NA
##
class Genome < GenomicData
  attr_reader :assembly
  
  # Initialize a 2bit file with a genomic reference sequence
  def initialize(filename)
    @data_file = File.expand_path(filename)
    @assembly = File.basename(@data_file, '.2bit')
    
    UCSCTools.twobit_info(@data_file).each { |chr,n| self[chr] = n }
  end
  
  # Check if this Genome includes the specified interval
  alias :super_include? :include?
  def include?(chr, start = nil, stop = nil)
    if start.nil?
      super_include?(chr)
    elsif stop.nil?
      super_include?(chr) and start >= 1 and start <= self[chr]
    else
      super_include?(chr) and start >= 1 and stop <= self[chr]
    end
  end
  
  # Get a specific stretch of sequence
  def sequence(chr, start = nil, stop = nil)
    raise GenomeError, "Genome does not include sequences for the range #{chr}:#{start}-#{stop}" if not include?(chr, start, stop)
    UCSCTools.twobit_to_fa(@data_file, chr, start, stop)
  end
  
  # Alias for sequence
  def query(chr, start = nil, stop = nil)
    sequence(chr, start, stop)
  end
  
  # Return the number of base pairs in the genome
  # Genome#length will return the number of chromosomes
  def num_bases
    self.map { |chr,n| n }.sum
  end
  
  # Reduce this genome to an Assembly object (just chromosome id's and their lengths)
  def to_assembly(name)
    a = Assembly.new(name, nil)
    self.each { |chr, n| a[chr] = n }
    return a
  end

  def to_s
    str = "Genome #{@assembly}: containing #{num_bases} base pairs\n"
    self.each do |chr, n|
      str += "\tChromosome #{chr} (length: #{n})\n"
    end
    
    return str
  end
end

class GenomeError < StandardError
end