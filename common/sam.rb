require 'entry_file'
require 'read_file'
require 'read'
require 'stringio'

# An entry in a SAM file
# From the specification, see: http://samtools.sourceforge.net/SAM-1.4.pdf
# Also helpful: http://chagall.med.cornell.edu/NGScourse/SAM.pdf
# TODO: Better handling of empty (*) values
class SAMEntry < Read
  attr_accessor :qname, :flag, :rname, :mapq, :cigar, :rnext, :pnext, :tlen
  
  def self.parse(line)
    record = line.chomp.split("\t")
  
    entry = self.new
    entry.qname = record[0]
    entry.flag = record[1].to_i
    entry.chr = record[2]
    entry.mapq = record[4].to_i
    entry.cigar = record[5]
    entry.rnext = record[6]
    entry.pnext = record[7].to_i
    entry.tlen = record[8].to_i
    entry.seq = Bio::Sequence::NA.new(record[9])
    entry.qual = record[10]
  
    extend = (entry.tlen == 0) ? entry.seq.length : entry.tlen.abs
  
    # According to the SAM specification, POS is the leftmost base
    # so 5' end for forward-mapping reads and
    # 3' end for reverse-complement mapping reads
    # Adjust start/stop appropriately
    if entry.watson?
      entry.start = record[3].to_i
      entry.stop = entry.start + extend - 1
    else
      entry.start = record[3].to_i + entry.seq.length - 1
      entry.stop = entry.start - extend + 1
    end

    return entry
  end
  
  def rname
    @chr
  end

  def pos
    @start
  end

  # FLAGS: See SAM specification
  def paired?
    (@flag & 0x0001) != 0
  end

  def unpaired?
    not paired?
  end

  def single?
    not paired?
  end

  def proper_pair?
    (@flag & 0x0002) != 0
  end

  def mapped?
    (@flag & 0x0004) == 0
  end

  def unmapped?
    not mapped?
  end

  def mate_mapped?
    (@flag & 0x0008) == 0
  end

  def mate_unmapped?
    not mate_mapped?
  end

  def watson?
    (@flag & 0x0010) == 0
  end

  def crick?
    not watson?
  end

  def mate_watson?
    (@flag & 0x0020) == 0
  end

  def mate_crick?
    not mate_watson?
  end

  def first?
    (@flag & 0x0040) != 0
  end

  def second?
    (@flag & 0x0080) != 0
  end

  def primary?
    (@flag & 0x0100) == 0
  end

  def failed?
    (@flag & 0x0200) != 0
  end

  def duplicate?
    (@flag & 0x0400) != 0
  end
end


# Access SAM files through stream operations
# Best option for iterating over large SAM files
class SAMFile < TextEntryFile
  extend ReadFile
  
  private
  
  def parse(line)
    SAMEntry.parse(line)
  end
end

# Access BAM files through stream operations using samtools
class BAMFile < BinaryEntryFile
  extend ReadFile
  
  def initialize(filename)
    super(filename)
    
    # Index the BAM file if its index can't be found
    index(@data_file) if not File.exist?(@data_file + '.bai')
  end
  
  # Cleanup the index
  def close
    File.delete(@data_file + '.bai') if File.exist?(@data_file + '.bai')
  end
  
  # Count the number of alignments in a given lookup
  def count(chr = nil, start = nil, stop = nil)
    %x[ samtools view -c #{@data_file} #{query_string(chr, start, stop)} ].chomp.to_i
  end

  private
  
  # Index sorted alignment for fast random access. Index file <aln.bam>.bai will be created.
  def index
    puts "Generating index for BAM file #{File.basename(@data_file)}" if ENV['DEBUG']
    %x[ samtools index #{@data_file} ]
  end
  
  def parse(line)
    SAMEntry.parse(line)
  end
  
  # Define how to query BAM files for lines
  def query_command(chr = nil, start = nil, stop = nil)
    "samtools view #{@data_file} #{query_string(chr, start, stop)}"
  end

  def self.query_string(chr = nil, start = nil, stop = nil)
    query = StringIO.new
    query << chr.to_s if chr
    query << ':' << start.to_s if start
    query << '-' << stop.to_s if stop

    return query.string
  end
end
