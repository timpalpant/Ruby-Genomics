require 'enumerator'
require 'unix_file_utils'
require 'chromosome'
require 'genomic_index_error'
require 'single_bp_data'
require 'single_bp_math'

##
# Lazy-load a Wig file, reading data into memory by chromosome
# only as needed to conserve memory
##
class WigFile
  include Enumerable
	include SingleBPMath
  
  attr_reader :name, :description, :data_file
  
  # Open a Wig file and parse its track/chromosome information
  def initialize(filename)
    @data_file = filename
    
    # Load the track information from the first line
    File.open(File.expand_path(filename)) do |f|
      track_line = f.gets.chomp
      track_line.split(' ').each do |opt|
        keypair = opt.split('=')
        key = keypair.first
        value = keypair.last[1..-2]
        
        # TODO: Load other track parameters
        case key
          when 'name'
            @name = value
          when 'description'
            @description = value
        end
      end
    end
    
    # Call grep to load the chromosome information:
    # Store hash of what chromosomes are available and what line they start at
    @index = Hash.new
		
		# Find fixedStep lines and index
		puts 'Indexing fixedStep lines' if ENV['DEBUG']
    File.grep_with_linenum(@data_file, 'fixedStep').each do |line|
      grep_line = line.split(':')
      line_num = grep_line.first.to_i
      fixedStep_line = grep_line.last
      
      fixedStep_line.split(' ').each do |opt|
        keypair = opt.split('=')
        key = keypair.first
        value = keypair.last
        
        if key == 'chrom'
          @index[value] = line_num
        end
      end
    end
		
		# Find variableStep lines and index
		puts 'Indexing variableStep lines' if ENV['DEBUG']
		File.grep_with_linenum(@data_file, 'variableStep').each do |line|
      grep_line = line.split(':')
      line_num = grep_line.first.to_i
      variableStep_line = grep_line.last
      
      variableStep_line.split(' ').each do |opt|
        keypair = opt.split('=')
        key = keypair.first
        value = keypair.last
        
        if key == 'chrom'
          @index[value] = line_num
        end
      end
    end
	
		# Raise an error if no chromosomes were found
		raise WigError, "No chromosome fixedStep/variableStep headers found in Wig file!" if @index.length == 0
  end
  
  # Enumerate over the chromosomes in this Wig file
  def each
    self.chromosomes.each do |chr_id| 
			yield [chr_id, self[chr_id]]
			
			# Hack to enforce garbage collection of old chromosomes
			# so that we don't exceed the available memory (hopefully)
			10.times { GC.start }
		end
		
		# Final GC (doing it multiple times seems to help?)
		10.times { GC.start }
  end

  # Return an array of all chromosomes in this WigFile file
  def chromosomes
    @index.keys
  end
	
	# Does this Wig file include data for chromosome chr?
	def include?(chr_id)
		@index.include?(chr_id)
	end
  
  # Allow indexing WigFiles like GenomicDatas
  def [](chr_id)
    chr(chr_id)
  end
    
  # Load data from disk and return a Vector of values for a given chromosome
  def chr(chr_id)
    # Call grep to find the chromosome within the Wig file
    raise GenomicIndexError, "Chromosome #{chr_id} not found in Wig file #{@data_file}!" unless include?(chr_id)
    start_line = @index[chr_id]
    
    # Read up to the next chromosome in the file
    end_line = @index.values.sort.select { |num| num > start_line }.first
    end_line -= 1 unless end_line.nil?
  
    # Call tail and head to get the appropriate lines from the Wig
    return Chromosome.load(@data_file, start_line, end_line)
  end
  
	# Output a summary about this Wig file
  def to_s
    str = "WigFile: connected to file #{@data_file}\n"
    @index.each do |chr,line|
      str += "\tChromosome #{chr} (starts at line: #{line})\n"
    end
    
    return str
  end
end