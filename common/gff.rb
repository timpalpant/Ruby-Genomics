require 'spot_array'
require 'spot'

##
# An entry 
##
class GFFEntry < Spot
	attr_accessor :chr
	
	def self.parse(line)
	  record = line.chomp.split("\t")
		
    spot = self.new
		spot.chr = record[0]
    spot.id = record[8]
    spot.start = record[3].to_i
    spot.stop = record[4].to_i
    spot.value = record[5].to_f
      
		return spot
	end
end

##
# Load GFF files completely into memory
##
class GFF < SpotArray
  # Load the GFF output file from NimbleScan
  def self.load(filename)
    puts "Loading GFF file: #{File.basename(filename)}" if ENV['DEBUG']
    spot_array = self.new
    
    File.foreach(filename) do |line|
      # Ignore comment lines
      next if line.start_with?('#') or line.chomp.empty?
      
      # Load spot lines
			entry = GFFEntry.parse(line)

			spot_array[entry.chr] ||= Array.new
      spot_array[entry.chr] << entry
    end
    
    puts "Loaded #{spot_array.num_spots} entries" if ENV['DEBUG']
    
    return spot_array
  end
end

##
# Stream GFF files by line or by chromosome
##
class GFFFile < File
	# Override each (line) to return each GFFEntry
	def self.foreach(filename)
		File.foreach(File.expand_path(filename)) do |line|
			# Skip comment and track lines
			next if line.start_with?('#') or line.start_with?('track')
			yield GFFEntry.parse(line)
		end
	end
end