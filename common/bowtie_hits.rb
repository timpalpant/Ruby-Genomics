require 'entry_file'
require 'read_file'
require 'read'

##
# Bowtie Output Files
# For the spec, see: http://bowtie-bio.sourceforge.net/manual.shtml#default-bowtie-output
##
class BowtieHitsEntry < Read
  attr_accessor :name, :refseq, :ceil, :mismatch
  
  def self.parse(line)
    entry = line.chomp.split("\t")
    
    read = self.new
    read.name = entry[0]
    strand = entry[1]
    read.chr = entry[2]
    read.start = entry[3].to_i
    read.seq = Bio::Sequence::NA.new(entry[4])
    read.stop = read.start + read.seq.length - 1
    read.qual = entry[5]
    read.ceil = entry[6]
    read.mismatch = entry[7]
    
    # Reverse-complement the read if it was on the - strand
    # Bowtie format specifies that the start reported is always the leftmost base
    # and the read sequence is reverse-complemented if it was on the - strand
    if strand == '-'
      tmp = read.start
      read.start = read.stop
      read.stop = tmp
      read.seq = read.seq.reverse_complement
      read.qual = read.qual.reverse
    end
    
    return read
  end
end

class BowtieHitsFile < TextEntryFile
  extend ReadFile

  private
  
  def parse(line)
    BowtieHitsEntry.parse(line)
  end
end