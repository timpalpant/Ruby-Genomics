require 'genomic_interval'
require 'genomic_data'

##
# Encapsulates information about an individual nucleosome
##
class Nucleosome < GenomicInterval
  attr_accessor :chr, :conditional_position, :dyad, :dyad_stdev, :dyad_mean, :occupancy

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
    nuke.occupancy = entry[7].to_i
    
    return nuke
  end
  
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
    
    NukeCallsFile.foreach(filename) do |nuke|
      calls[nuke.chr] ||= Array.nwe
      calls[nuke.chr] << nuke
    end
    
    return calls
  end
end

##
# Lists of Nucleosome Calls
##
class NukeCallsFile
  # Override each (line) to return each Nucleosome entry
	def self.foreach(filename, chr = nil)
    if chr.nil?
      File.foreach(File.expand_path(filename)) do |line|
        # Skip comment lines
        next if line.start_with?('#') or line.chomp.empty?
        yield Nucleosome.parse(line)
      end
    else
      IO.popen("grep -w #{chr} #{File.expand_path(filename)}") do |output|
        output.each do |line|
          # Skip comment lines
          next if line.start_with?('#') or line.chomp.empty?
          yield Nucleosome.parse(line)
        end
      end
    end
	end
end