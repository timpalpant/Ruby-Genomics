require 'entry_file'
require 'spot_file'
require 'spot'

##
# An entry in a Bed file
##
class BedEntry < Spot 
  attr_accessor :thick_start, :thick_end, :item_rgb, :block_count, :block_sizes, :block_starts

  def self.parse(line)
    begin
      entry = line.chomp.split("\t")
      raise BedError, "Invalid Bed Entry: Bed must have at least 3 columns" if entry.length < 3
        
      spot = self.new
      spot.chr = entry[0]
      # Bed is 0-indexed
      spot.start = entry[1].to_i+1
      # And half open
      spot.stop = entry[2].to_i  
      spot.id = entry[3] if entry.length >= 4
      spot.value = entry[4].to_f if entry.length >= 5
      
      # Reverse start/stop if on the - strand
      if entry.length >= 6 and entry[5].chomp == '-' and spot.start < spot.stop
        tmp = spot.start
        spot.start = spot.stop
        spot.stop = tmp
      end
  
      if entry.length >= 8
        spot.thick_start = entry[6].to_i + 1
        spot.thick_end = entry[7].to_i
      end

      spot.item_rgb = entry[8] if entry.length >= 9

      if entry.length >= 12
        spot.block_count = entry[9].to_i 
        spot.block_sizes = entry[10].split(',').map { |v| v.to_i }
        raise BedError, "Invalid Bed Entry: blockCount does not correspond to number of blockSizes" if spot.block_count != spot.block_sizes.length
        spot.block_starts = entry[11].split(',').map { |v| spot.low + v.to_i }
        raise BedError, "Invalid Bed Entry: blockCount does not correspond to number of blockStarts" if spot.block_count != spot.block_starts.length
      end

      return spot
    rescue
      raise BedError, "Invalid Bed Entry!"
    end
  end
  
  def name
    @id
  end

  def to_bed6
    "#{@chr}\t#{low-1}\t#{high}\t#{name ? name : '.'}\t#{@value ? @value : '.'}\t#{strand}"
  end

  def to_bed12
    s = to_bed6
    s += "\t#{@thick_start ? @thick_start-1 : 0}\t#{@thick_end ? @thick_end : 0}\t#{@item_rgb ? @item_rgb : 0}"
    s += "\t#{@block_count}\t#{@block_sizes.join(',')}\t#{@block_starts.map { |v| v-low }.join(',')}" if @block_count
    return s
  end

  def to_bed
    # Write Bed-12 fields if they are defined
    if @thick_start or @thick_end or @item_rgb or @block_count or @block_sizes or @block_starts
      to_bed12
    else
      to_bed6
    end
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
