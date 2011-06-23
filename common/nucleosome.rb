require 'entry_file'
require 'spot_file'
require 'spot'

##
# Encapsulates information about an individual nucleosome
##
class Nucleosome < Spot
  attr_accessor :conditional_position, :dyad, :dyad_stdev, :dyad_mean

  # Parse a NukeCallsFile entry
  def self.parse(line)
    begin
      entry = line.chomp.split("\t")
      
      nuke = Nucleosome.new
      nuke.chr = entry[0]
      nuke.start = entry[1].to_i
      nuke.stop = entry[2].to_i
      nuke.dyad = entry[3].to_i
      nuke.dyad_stdev = entry[4].to_f
      nuke.conditional_position = entry[5].to_f
      nuke.dyad_mean = entry[6].to_i
      nuke.value = entry[7].to_f
      
      return nuke
    rescue
      raise NucleosomeError, "Invalid nucleosome call entry!"
    end
  end

  def occupancy=(value)
    @value = value
  end
  
  def occupancy
    @value
  end
  
  # Use the dyad as the nucleosome position
  def position
    dyad
  end
  
  def to_s
    "#{@chr}\t#{@start}\t#{@stop}\t#{@dyad}\t#{@dyad_stdev}\t#{@conditional_position}\t#{@dyad_mean}\t#{@occupancy}"
  end
end

##
# Lists of Nucleosome Calls
##
class NukeCalls
  HEADER = "#Chromosome\tNuke Start\tNuke Stop\tDyad\tDyad StDev\tSmoothed Position\tDyad Mean\tDyad Count"
end

##
# Lists of Nucleosome Calls
##
class NukeCallsFile < TextEntryFile
  include SpotFile
 
  CHR_COL = 1
  START_COL = 4
  END_COL = 4
  
  def initialize(filename)
    super(filename, CHR_COL, START_COL, END_COL)
  end

  private
  
  def parse(line)
    Nucleosome.parse(line)
  end
end

class NucleosomeError < EntryFileError
end
