require 'entry_file'
require 'spot_file'
require 'spot'

##
# An entry in a bedGraph file
##
class BedGraphEntry < Spot
  def self.parse(line)
    begin
      entry = line.chomp.split("\t")
        
      spot = self.new
      spot.chr = entry[0]
      spot.start = entry[1].to_i
      spot.stop = entry[2].to_i
      spot.value = entry[3].to_f if entry.length >= 4
        
      return spot
    rescue
      raise BedGraphError, "Invalid BedGraph Entry!"
    end
  end
end

##
# Stream bedgraph files by line or by chromosome
##
class BedGraphFile < TextEntryFile
  include SpotFile
  
  CHR_COL = 1
  START_COL = 2
  END_COL = 3
  
  def initialize(filename)
    super(filename, CHR_COL, START_COL, END_COL)
  end

  private
  
  # Define how to parse BedGraph entries
  def parse(line)
    BedGraphEntry.parse(line)
  end
end

class BedGraphError < EntryFileError
end