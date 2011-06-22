require 'entry_file'
require 'spot_file'
require 'spot'

##
# An entry in a bedGraph file
##
class BedGraphEntry < Spot
  def self.parse(line)
    entry = line.chomp.split("\t")
      
    raise BedGraphError, "Invalid BedGraph entry!" if entry.length < 3
      
    spot = self.new
    spot.chr = entry[0]
    spot.start = entry[1].to_i
    spot.stop = entry[2].to_i
    spot.value = entry[3].to_f if entry.length >= 4
      
    return spot
  end
end

##
# Stream bedgraph files by line or by chromosome
##
class BedGraphFile < TextEntryFile
  extend SpotFile

  # Convert a bedGraph file to a BigWig
  def self.to_bigwig(input_file, output_file, assembly)
    # BedGraph must be sorted first
    tmp_sorted = File.expand_path(input_file + '.sorted')
    File.sort(File.expand_path(input_file), tmp_sorted, '-k1,1 -k2,2')
    %x[ bedGraphToBigWig #{tmp_sorted} #{File.expand_path(assembly.len_file)} #{File.expand_path(output_file)} ]
    File.delete(tmp_sorted)
  end

  private
  
  # Define how to parse BedGraph entries
  def parse(line)
    BedGraphEntry.parse(line)
  end
end

class BedGraphError < StandardError
end