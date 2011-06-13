require 'spot_array'
require 'spot'

##
# An entry in a MACS output file
##
class MACSEntry < Spot
	attr_accessor :chr, :summit, :tags, :fold_enrichment, :fdr
	
	def self.parse(line)
    entry = line.chomp.split("\t")
			
    spot = self.new
		spot.chr = entry[0]
    spot.start = entry[1].to_i
    spot.stop = entry[2].to_i
    spot.summit = spot.start + entry[4].to_i
    spot.tags = entry[5].to_i
    spot.fold_enrichment = entry[7].to_f
      
    return spot
	end
end

##
# Load MACS output files
##
class MACS < SpotArray
  # Load the Bed file as a SpotArray
  def self.load(filename)
    puts "Loading MACS file: #{File.basename(filename)}" if ENV['DEBUG']
    spot_array = self.new
    
    File.foreach(File.expand_path(filename)) do |line|
      # Ignore comment lines
      next if line.start_with?('#') or line.chomp.empty?
      
			# Parse the entry
			entry = MACSEntry.parse(line)
			
      spot_array[entry.chr] ||= Array.new
      spot_array[entry.chr] << entry
    end

    puts "Loaded #{spot_array.num_spots} peak calls" if ENV['DEBUG']
    return spot_array
  end
end
