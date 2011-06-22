require 'entry_file'
require 'spot_file'

##
# Encapsulates information about an individual nucleosome
##
class Nucleosome < Spot
  attr_accessor :conditional_position, :dyad, :dyad_stdev, :dyad_mean

  # Parse a NukeCallsFile entry
  def self.parse(line)
    entry = line.chomp.split("\t")
    raise "Not a valid nucleosome call entry!" if entry.length < 8
    
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
  extend SpotFile

  private
  
  def parse(line)
    Nucleosome.parse(line)
  end
end
