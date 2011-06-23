require 'entry_file'
require 'spot_file'
require 'spot'

##
# An entry in a Bed file
##
class BedEntry < Spot 
  attr_accessor :strand, :thick_start, :thick_end, :item_rgb, :block_count, :block_sizes, :block_starts

  def self.parse(line)
    begin
      entry = line.chomp.split("\t")
        
      spot = self.new
      spot.chr = entry[0]
      spot.start = entry[1].to_i
      spot.stop = entry[2].to_i
      # Bed format specifies that start must be less than or equal to stop
      raise BedError, "Invalid Bed entry: start > stop" if spot.start > spot.stop
      spot.id = entry[3] if entry.length >= 4
      spot.value = entry[4].to_f if entry.length >= 5
      
      # Reverse start/stop if on the - strand
      strand = entry[5] if entry.length >= 6
      if strand == '-'
        tmp = start.start
        spot.start = spot.stop
        spot.stop = tmp
      end
  
      spot.thick_start = entry[6] if entry.length >= 7
      spot.thick_end = entry[7] if entry.length >= 8
      spot.item_rgb = entry[8] if entry.length >= 9
      spot.block_count = entry[9] if entry.length >= 10
      spot.block_sizes = entry[10] if entry.length >= 11
      spot.block_starts = entry[11] if entry.length >= 12
      
      return spot
    rescue
      raise BedError, "Invalid Bed Entry!"
    end
  end
  
  def name
    @id
  end

  def to_bed    
    # Write Bed-6 format by default
    s = "#{@chr}\t#{low}\t#{high}\t#{name}\t#{@value}\t#{strand}"
    
    # Write Bed-12 fields if they are defined
    if @thick_start or @thick_end or @item_rgb or @block_count or @block_sizes or @block_starts
      s += "\t#{@thick_start}\t#{@thick_end}\t#{@item_rgb}\t#{@block_count}\t#{@block_sizes}\t#{@block_starts}"
    end
    
    return s
  end
end

##
# Get data from BigBed files
##
class BigBedFile < BinaryEntryFile
  def initialize(filename)
    raise BedError, "BigBed files are not yet implemented"
  end
end

##
# Stream bed files by line or by chromosome
##
class BedFile < TextEntryFile
  include SpotFile
  
  CHR_COL = 1
  START_COL = 2
  END_COL = 3
  
  def initialize(filename)
    super(filename, CHR_COL, START_COL, END_COL)
  end

  private
  
  # Define how to parse Bed entries
  def parse(line)
    BedEntry.parse(line)
  end
end

class BedError < EntryFileError
end
