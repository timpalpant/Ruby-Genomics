require 'spot_array'
require 'spot'

##
# An entry in a Bed file
##
class BedEntry < Spot
	attr_accessor :chr
	
	def self.parse(line)
		entry = line.chomp.split("\t")
      
		raise BedError, "Invalid Bed Entry!" if entry.length < 3
			
    spot = self.new
		spot.chr = entry[0]
    spot.start = entry[1].to_i
    spot.stop = entry[2].to_i
    spot.id = entry[3] if entry.length >= 4
    spot.value = entry[4].to_f if entry.length >= 5
		
		return spot
	end
end

##
# Load Bed files completely into memory
# @DEPRECATED: only feasible for small genomess / datasets
##
class Bed < SpotArray
  # Load the Bed file as a SpotArray
  def self.load(filename)
		puts "Loading Bed file: #{File.basename(filename)}" if ENV['DEBUG']
    spot_array = self.new
		
    BedFile.foreach(filename) do |entry|
      spot_array[entry.chr] ||= Array.new
      spot_array[entry.chr] << entry
    end

    puts "Loaded #{spot_array.num_spots} entries" if ENV['DEBUG']
    return spot_array
  end
end

##
# Base class for all Bed files (BedFile / BigBedFile)
##
class AbstractBedFile
end

##
# Get data from BigBed files
##
class BigBedFile < AbstractBedFile
end

##
# Stream bed files by line or by chromosome
##
class BedFile < AbstractBedFile
	# Override each (line) to return each BedEntry
	def self.foreach(filename)
		File.foreach(File.expand_path(filename)) do |line|
			# Skip comment and track lines
			next if line.start_with?('#') or line.start_with?('track') or line.chomp.empty?
			yield BedEntry.parse(line)
		end
	end
end

class BedError < StandardError
end
