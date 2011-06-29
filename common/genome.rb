require 'bio'
require 'assembly'
require 'ucsc_tools'
require 'stats'
require 'enumerator'
require 'stringio'

##
# Get genomic sequences from 2bit files
# Loads the sequences as Bio::Sequence::NA
##
class Genome < GenomicData
  # Initialize a 2bit file with a genomic reference sequence
  def initialize(filename)
    @data_file = File.expand_path(filename)
    @assembly = File.basename(@data_file, '.2bit')
    
    UCSCTools.twobit_info(filename).each do |chr,n|
      self[chr] = n
    end
  end
  
  # Iterate over the chromosomes with their lengths like an Assembly
  def each
    @index.each do |chr,n|
      yield(chr,n)
    end
  end
  
  # Get a specific stretch of sequence
  def sequence(chr, start = nil, stop = nil)
    s = StringIO.new
    UCSCTools.twobit_to_fa(@data_file, chr, start, stop) do |line|
      # Skip the header line
      next if line.start_with?('>')
      
      s << line.chomp
    end
    
    Bio::Sequence::NA.new(s.string)
  end
  
  # Alias for sequence
  def query(chr, start = nil, stop = nil)
    sequence(chr, start, stop)
  end
  
  # Return the number of base pairs in the genome
  # Genome#length will return the number of chromosomes
  def num_bases
    @index.map { |chr,n| n }.sum
  end
  
  # Reduce this genome to an Assembly object (just chromosome id's and their lengths)
  def to_assembly(name)
    a = Assembly.new(name, nil)
    
    self.each do |chr, n|
      a[chr] = n
    end
    
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
