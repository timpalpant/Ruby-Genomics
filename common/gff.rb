require 'entry_file'
require 'spot_file'
require 'spot'

##
# An entry in a GFF file
# For the format spec, see: http://genome.ucsc.edu/FAQ/FAQformat.html#format3
##
class GFFEntry < Spot
  attr_accessor :source, :feature, :frame
  
  def self.parse(line)
    begin
      record = line.chomp.split("\t")
      
      spot = self.new
      spot.chr = record[0]
      spot.source = record[1]
      spot.feature = record[2]
      spot.start = record[3].to_i
      spot.stop = record[4].to_i
      spot.value = record[5].to_f
      strand = record[6]
      spot.frame = record[7]
      spot.id = record[8]
      
      # Ensure that the coordinates (start/stop) match the strand, if specified
      if strand == '+'
        tmp_low = spot.low
        tmp_high = spot.high
        spot.start = tmp_low
        spot.stop = tmp_high
      elsif strand == '-'
        tmp_low = spot.low
        tmp_high = spot.high
        spot.start = tmp_high
        spot.stop = tmp_low
      end
        
      return spot
    rescue
      raise GFFError, "Not a valid GFF Entry" 
    end
  end
  
  def seqname
    @chr
  end
  
  def group
    @id
  end
  
  def score
    @value
  end
end


##
# Stream GFF files
##
class GFFFile < TextEntryFile
  extend SpotFile
  
  CHR_COL = 1
  START_COL = 4
  END_COL = 5
  
  def initialize(filename)
    super(filename, CHR_COL, START_COL, END_COL)
  end

  private
  
  # Define how to parse GFF entries
  def parse(line)
    GFFEntry.parse(line)
  end
end


class GFFError < EntryFileError
end