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
##
class Bed < SpotArray
  # Load the Bed file as a SpotArray
  def self.load(filename)
		puts "Loading Bed file: #{File.basename(filename)}" if ENV['DEBUG']
    spot_array = self.new
    skipped = 0
		
    File.foreach(File.expand_path(filename)) do |line|
      # Ignore comment/empty lines
      next if line.start_with?('#') or line.start_with?('track') or line.chomp.empty?
      
      # Load spot lines
			begin
				entry = BedEntry.parse(line)
			rescue BedError
				puts "Invalid line: #{line}"
				skipped += 1
				next
			end
			
      spot_array[entry.chr] ||= Array.new
      spot_array[entry.chr] << entry
    end

    puts "Loaded #{spot_array.num_spots} entries" if ENV['DEBUG']
		puts "Skipped #{skipped} invalid entries" if ENV['DEBUG'] and skipped > 0
    
    return spot_array
  end
end

##
# Stream bed files by line or by chromosome
##
class BedFile < File
	# Override each (line) to return each BedEntry
	def self.foreach(filename)
		File.foreach(File.expand_path(filename)) do |line|
			# Skip comment and track lines
			next if line.start_with?('#') or line.start_with?('track')
			yield BedEntry.parse(line)
		end
	end
end

class BedError < StandardError
end