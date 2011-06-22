require 'entry_file'
require 'spot_file'
require 'spot'

##
# An entry in a Bed file
##
class BedEntry < Spot 
  attr_accessor :strand, :thick_start, :thick_end, :item_rgb, :block_count, :block_sizes, :block_starts

  def self.parse(line)
    entry = line.chomp.split("\t")
      
    raise BedError, "Invalid Bed Entry!" if entry.length < 3
      
    spot = self.new
    spot.chr = entry[0]
    spot.start = entry[1].to_i
    spot.stop = entry[2].to_i
    spot.id = entry[3] if entry.length >= 4
    spot.value = entry[4].to_f if entry.length >= 5
    spot.strand = entry[5] if entry.length >= 6
    spot.thick_start = entry[6] if entry.length >= 7
    spot.thick_end = entry[7] if entry.length >= 8
    spot.item_rgb = entry[8] if entry.length >= 9
    spot.block_count = entry[9] if entry.length >= 10
    spot.block_sizes = entry[10] if entry.length >= 11
    spot.block_starts = entry[11] if entry.length >= 12
    
    return spot
  end
  
  def name
    @id
  end
end

##
# Get data from BigBed files
##
class BigBedFile < BinaryEntryFile
end

##
# Stream bed files by line or by chromosome
##
class BedFile < TextEntryFile
  extend SpotFile

  private
  
  # Define how to parse Bed entries
  def parse(line)
    BedEntry.parse(line)
  end
end

class BedError < StandardError
end
