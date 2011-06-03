require 'spot_array'
require 'spot'

##
# An entry in a bedGraph file
##
class BedGraphEntry < Spot
	attr_accessor :chr

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
# Load bedGraph files
# @DEPRECATED: only feasible for small genomess / datasets
##
class BedGraph < SpotArray
  # Load the Bed file as a SpotArray
  def self.load(filename)
		puts "Loading BedGraph file: #{File.basename(filename)}" if ENV['DEBUG']
    spot_array = self.new
    
		skipped = 0
    File.foreach(File.expand_path(filename)) do |line|
      # Ignore comment lines
      next if line[0] == '#' or line.start_with?('track') or line.chomp.empty?
      
      # Load spot lines
			begin
				entry = BedGraphEntry.parse(line)
			rescue BedGraphError
				skipped += 1
				puts "Skipped line: #{line}" if ENV['DEBUG']
				next
			end
			
			spot_array[entry.chr] ||= Array.new
      spot_array[entry.chr] << entry
    end

    puts "Loaded #{spot_array.num_spots} entries" if ENV['DEBUG']
		puts "Skipped #{skipped} invalid lines" if ENV['DEBUG'] and skipped > 0
    
    return spot_array
  end
  
  # Convert a bedGraph file to a BigWig
  def self.to_bigwig(input_file, output_file, assembly)
    # BedGraph must be sorted first
    tmp_sorted = input_file + '.sorted'
    File.sort(input_file, tmp_sorted, '-k1,1 -k2,2')
    %x[ bedGraphToBigWig #{tmp_sorted} #{assembly.len_file} #{output_file} ]
    File.delete(tmp_sorted)
  end
end

##
# Stream bedgraph files by line or by chromosome
##
class BedGraphFile < File
	# Override each (line) to return each BedGraphEntry
	def self.foreach(filename)
		File.foreach(File.expand_path(filename)) do |line|
			# Skip comment and track lines
			next if line.start_with?('#') or line.start_with?('track')
			yield BedGraphEntry.parse(line)
		end
	end
end

class BedGraphError < StandardError
end