require 'genomic_interval'
require 'genomic_data'

##
# Encapsulates information about an individual nucleosome
##
class Nucleosome < GenomicInterval
  attr_accessor :conditional_position, :dyad, :dyad_stdev, :dyad_mean, :occupancy

  # Use the dyad as the nucleosome position
  def position
    dyad
  end
  
  def to_s
    "#{@start}\t#{@stop}\t#{@dyad}\t#{@dyad_stdev}\t#{@conditional_position}\t#{@dyad_mean}\t#{@occupancy}"
  end
end

##
# Lists of Nucleosome Calls
##
class NukeCalls < GenomicData
	
	HEADER = "#Chromosome\tNuke Start\tNuke Stop\tDyad\tDyad StDev\tSmoothed Position\tDyad Mean\tDyad Count"
	
  # Load a list of nucleosome calls
  def self.load(filename)
    calls = self.new
    
    File.foreach(filename) do |line|
      # Skip header lines
      next if line[0] == '#'
      entry = line.split("\t")
      
      nuke = Nucleosome.new
      nuke.start = entry[1].to_i
      nuke.stop = entry[2].to_i
      nuke.dyad = entry[3].to_i
      nuke.dyad_stdev = entry[4].to_f
      nuke.conditional_position = entry[5].to_f
      nuke.dyad_mean = entry[6].to_i
      nuke.occupancy = entry[7].to_i
      
			chr = entry[0]
			calls[chr] ||= Array.new
      calls[chr] << nuke
    end
    
    return calls
  end
end