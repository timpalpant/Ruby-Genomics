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
      break unless track_line.start_with?('track')
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
		
		# Find chromosome header lines and index (ghetto B-tree)
		puts 'Indexing chromosome header lines' if ENV['DEBUG']
    File.grep_with_linenum(@data_file, 'chrom').each do |line|
      grep_line = line.split(':')
      line_num = grep_line.first.to_i
      header_line = grep_line.last
      next unless header_line.start_with?('fixedStep') or header_line.start_with?('variableStep')
      
      header_line.split(' ').each do |opt|
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
			yield [chr_id, chr(chr_id)]
		end
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
    raise GenomicIndexError, "Chromosome #{chr_id} not found in Wig file #{@data_file}!" unless include?(chr_id)
  
    # Call tail and head to get the appropriate lines from the Wig
    return Chromosome.load(@data_file, chr_start(chr_id), chr_stop(chr_id))
  end
  
  # Get the starting line for a chromosome
  def chr_start(chr_id)
    raise GenomicIndexError, "Chromosome #{chr_id} not found in Wig file #{@data_file}!" unless include?(chr_id)
    @index[chr_id]
  end
  
  # Get the stop line for a chromosome
  def chr_stop(chr_id)
    raise GenomicIndexError, "Chromosome #{chr_id} not found in Wig file #{@data_file}!" unless include?(chr_id)
    start_line = chr_start(chr_id)
    
    # Read up to the next chromosome in the file, or the end if there are no more chromosomes
    end_line = @index.values.sort.select { |num| num > start_line }.first
    end_line -= 1 unless end_line.nil?
    end_line = File.num_lines(@data_file) if end_line.nil?
    
    return end_line
  end
  
  # Get the length of a chromosome from the index
  def chr_length(chr_id)
    raise GenomicIndexError, "Chromosome #{chr_id} not found in Wig file #{@data_file}!" unless include?(chr_id)
    chr_stop(chr_id) - chr_start(chr_id)
  end
  
  # Return single-bp data from the specified region
  def query(chr, start, stop)
    raise GenomicIndexError, "Chromosome #{chr} not found in Wig file #{@data_file}!" unless include?(chr)
    
    # Parse the header of the desired chromosome
    header = File.lines(@data_file, chr_start(chr), chr_start(chr)).first
    raise GenomicIndexError, 'Random queries are only available for fixedStep-style chromosomes' unless header.start_with?('fixedStep')
    parsed = Chromosome.parse_header(header)
    
    raise 'Random queries are not yet implemented for data with step != 1' if parsed.step != 1
    raise GenomicIndexError, 'Specified interval outside of data range in Wig file!' if start < parsed.start or stop > parsed.start + parsed.step*chr_length(chr)

    # Calculate the lines needed
    low = [start, stop].min
    high = [start, stop].max
    start_line = chr_start(chr) + 1 + (low - parsed.start)/parsed.step
    stop_line = start_line + (high-low)/parsed.step
    
    # Call tail and head to get the appropriate lines from the Wig
    values = File.lines(@data_file, start_line, stop_line).map { |line| line.to_f }
    values.reverse! if start > stop
    
    return values
  end
  
	# Output a summary about this Wig file
  def to_s
    str = "WigFile: connected to file #{@data_file}\n"
    @index.each do |chr,line|
      str += "\tChromosome #{chr} (lines: #{line}..#{chr_stop(chr)})\n"
    end
    
    return str
  end
end
