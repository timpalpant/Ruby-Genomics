require 'spot_array'

##
# Load Affymetrix array files
# @DEPRECATED: only feasible for small genomess / datasets
##
class Affy < SpotArray
  # Load the Affy file as a SpotArray
  def self.load(filename)
		puts "Loading Affymetrix file: #{File.basename(filename)}" if ENV['DEBUG']
    spot_array = self.new
    
		prev_chr, prev_start, prev_value = nil, nil, nil
		File.foreach(File.expand_path(filename)) do |line|
			# Skip comment lines
			next if line.start_with?('#')
    
			entry = line.chomp.split("\t")
			chr = entry[0]
			start = entry[1].to_i
			value = entry[2].to_f
    
			# Extend values (spots) up to the next spot (i.e. tile)
			if prev_chr and prev_start and prev_value and prev_chr == chr and prev_start < start
				spot = Spot.new
				spot.start = prev_start
				spot.stop = start-1
				spot.value = prev_value
				spot_array[chr] ||= Array.new
				spot_array[chr] << spot
			end
		
			prev_chr = chr
			prev_start = start
			prev_value = value
		end

    puts "Loaded #{spot_array.num_spots} spots" if ENV['DEBUG']
    
    return spot_array
  end
end