require 'spot'
require 'stats'
require 'genomic_data'
require 'assembly'
require 'wig'

##
# A GenomicData with values for each Spot, i.e. a micrarray dataset
##
class SpotArray < GenomicData
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
  
	# Override Array#uniq, taking the median value of replicate spots
	# TODO: override Array#uniq!
	def uniq
		uniq_spots = self.inject({}) do |hash, spot| 
		  hash[spot] ||= []
		  hash[spot] << spot.value
		  hash
		end
		
		uniq_spots.inject(SpotArray.new) do |spots, spot| 
			spot.first.value = spot.last.median
		  spots << spot.first
		end
	end
	
	# Write this array to GFF
	def to_gff(filename)
    File.open(filename,'w') do |f|
      self.each do |chr,spots|
      	spots.each do |spot|
      		f.puts "{chr}\t" + spot.to_gff
      	end
      end
    end
	end
	
	#  Write this array to bed format
	def to_bed(filename)
    File.open(filename,'w') do |f|
  		self.each do |chr,spots|
  			spots.each do |spot|
  				f.puts "{chr}\t" + spot.to_bed
  			end
  		end
  	end
	end
	
	#  Write this array to bedGraph format
	def to_bedGraph(filename)
    File.open(filename,'w') do |f|
  		self.each do |chr,spots|
  			spots.each do |spot|
  				f.puts "{chr}\t" + spot.to_bedGraph
  			end
  		end
  	end
	end
	
	# Write this array to Wig format
	# Construct a Wig that is the data from all spots
	# averaged, if multiple spots cover a given base
	def to_wig(filename, assembly = Assembly.yeast)		
		# Iterate over each chromosome, mapping all spots and averaging
		File.open(filename,'w') do |f|
			# TODO: should be rewritten to intelligently use step size
      f.puts Wig.track_header(@name, @description) 
			
			self.each do |chr,spots|
				# Skip if this chromosome is not in the specified assembly
				next unless assembly.include?(chr)
				
				# Allocate space for the new Wig chromosomes
				values = Chromosome.new(assembly[chr])
				count = Chromosome.new(assembly[chr])
				
				spots.each do |spot|
					# Get the high and low spot coordinates, and clamp to the ends of the chromosome
					low = [1, spot.low].max
					high = [spot.high, values.length].min
					
					# Skip if the spot is completely outside of the chromosome
					next if low >= high
					
					# Add the spot data to the wig
					for bp in low-1..high-1
                                          values[bp] += spot.value
                                          count[bp] += 1
                                        end
				end
				
				# Average spots that overlap (NaN where there isn't any data)
        for i in 0...values.length
          values[i] /= count[i] unless count[i] == 0
        end
			
				# Write to output file
				f.puts Wig.fixed_step(chr, values)
        f.puts values
			end
    end
	end
				
	# Return a value for the given location
	def value(chr, base)
		spots = self[chr].select { |spot| spot.include?(base) }
			
		# If only one spot covers the given location, return its value
		if spots.length == 1
			return spots.first.value
		# If multiple spots cover the location, return their average
		elsif spots.length > 1
			return spots.collect { |spot| spot.value }.mean
		# Otherwise return nil
		else
			return nil
		end
	end
end
