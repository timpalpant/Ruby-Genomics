require 'spot_array'
require 'roman_numerals'

# Load the MasterAnnotations data into a Hash that maps
# Systematic ID (ORF or Intergenic) -> Data
# Similar to VLOOKUP in Excel
class MasterAnnotations < Hash

	# Create a new master annotations database
	def initialize
		f = File.open("../resources/MasterAnnotations.txt")
		header = f.gets.split("\t")
		
		while(line = f.gets)
			record = line.split("\t")
			name = record.first
			self[name] = record[1..-1]
		end
	end
	
	# Load master annotations database as a SpotArray
	# with sacCer2 coordinates
	def self.to_spot_array
		annotations = SpotArray.new
		
		f = File.open("../resources/MasterAnnotations.txt")
		header = f.gets.split("\t")
		
		while(line = f.gets)
			record = line.split("\t")
			
			spot = Spot.new
			spot.id = record[0]
			chr = record[9].to_i.to_roman
			spot.start = record[10].to_i
			spot.stop = record[11].to_i
			
			annotations[chr] ||= Array.new
			annotations[chr] << spot
		end
		
		return annotations
	end
end