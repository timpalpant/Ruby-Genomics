require 'spot'
require 'stats'
require 'genomic_data'
require 'assembly'
require 'wig'

##
# A GenomicData with values for each Spot, i.e. a micrarray dataset
##
class SpotArray < GenomicData

  ##
  # STATISTICAL METHODS
  ##

	# The total number of spots in this SpotArray
  def num_spots
    self.collect { |chr_id,spots| spots.length }.sum
  end
	
	# The sum of the values of all spots
	def total
		self.collect { |chr_id,spots| spots.collect { |spot| spot.value }.sum }.sum
	end
	
	# The mean value of all spots
	def mean
		total.to_f / num_spots
	end
	
	# The standard deviation of all spots
	def stdev(mean = self.mean)
		self.map { |chr_id,spots| spots.map { |spot| (spot.value-mean)**2 }.sum }.sum / num_spots
	end
	
	##
	# QUERY METHODS
	##
	
	# Return values for the given window, with single-bp resolution (even if inferred)
	def query(chr, start, stop)
	  low = [start, stop].min
	  high = [start, stop].max
	  length = high - low + 1
	  
	  total = Array.new(length, 0)
	  count = Array.new(length, 0)
	  
	  self[chr].select { |spot| spot.high >= low and spot.low <= high }.each do |spot|
			# Get the high and low spot coordinates, and clamp to the ends of the window
			low = [low, spot.low].max
			high = [spot.high, high].min
	  
	    for bp in spot.low..spot.high
	      total[bp-low] += spot.value unless spot.value.nil?
	      count[bp-low] += 1
	    end
	  end
	  
	  # Map base pairs without data to nil, and take the mean of overlapping probes
	  med = Array.new(length) do |i|
	    if count[i] > 0
	      total[i].to_f / count[i]
	    else
	      nil
	    end
	  end
	  
	  # Allow Crick querying
	  med.reverse! if start > stop
	  return med
	end
				
	# Return a value for the given location
	def value(chr, base)
		spots = self[chr].select { |spot| spot.include?(base) }
			
		# If only one spot covers the given location, return its value
		if spots.length == 1
			return spots.first.value
		# If multiple spots cover the location, return their median
		elsif spots.length > 1
			return spots.collect { |spot| spot.value }.median
		# Otherwise return nil
		else
			return nil
		end
	end
  
  ##
  # OUTPUT-TO-DISK METHODS
  ##
	
	# Write this array to GFF
	def to_gff(filename)
    self.to_disk(filename) do |chr,spot|
      "#{chr}\t" + spot.to_gff
    end
	end
	
	#  Write this array to bed format
	def to_bed(filename)
    self.to_disk(filename) do |chr,spot|
      "#{chr}\t" + spot.to_bed
    end
	end
	
	#  Write this array to bedGraph format
	def to_bedGraph(filename)
    self.to_disk(filename) do |chr,spot|
      "#{chr}\t" + spot.to_bedGraph
    end
	end
	
	# Write this array to Wig format
	# Construct a Wig that is the data from all spots
	# averaged, if multiple spots cover a given base
	def to_wig(filename, assembly)		
		# Iterate over each chromosome, mapping all spots and averaging
		File.open(filename,'w') do |f|
			# TODO: should be rewritten to intelligently use step size
      f.puts Wig.track_header(@name, @description) 
			
			self.chromosomes.each do |chr|
				# Skip if this chromosome is not in the specified assembly
				next unless assembly.include?(chr)
				
				# Allocate space for the new Wig chromosomes
				values = query(chr, 1, assembly[chr]).to_chr(1, 1, 1)
			
				# Write to output file
				f.puts Wig.fixed_step(chr, values)
        f.puts values
			end
    end
	end
	
	private
	
	def to_disk(filename)
	  File.open(File.expand_path(filename), 'w') do |f|
      self.each do |chr,spots|
      	spots.each do |spot|
      		f.puts yield(chr, spot)
      	end
      end
    end
	end
end
