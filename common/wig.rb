require 'enumerator'
require 'unix_file_utils'
require 'chromosome'
require 'genomic_index_error'
require 'stats'
require 'genomic_data'
require 'assembly'
require 'gsl'

##
# Load a fixedStep Wig file completely into memory
# @DEPRECATED: only feasible for small genomes / datasets
##
class Wig	
  attr_accessor :name, :description, :assembly
  
  def self.load(filename)
    wig = self.new
    
		# Piggyback on WigFile parsing
    # Iterate over each chromosome and store permanently
    cached_wig = WigFile.new(filename)
    wig.name = cached_wig.name
    wig.description = cached_wig.description
    
    cached_wig.each do |chr,values|
    	wig[chr] = values
    end
    
    return wig
  end
  
  # Return a new empty SingleBPData (zeros) with dimensions of an assembly
  def self.for_assembly(a)    
    empty_seq = self.new
		empty_seq.assembly = a.name
  	a.each { |chr,num_bases| empty_seq[chr] = Chromosome.new(num_bases) }
    return empty_seq
  end
  
  # Summary information about this sequencing file object
  def to_s
    out = "Sequencing File: #{@name}\n"
    out << "Assembly: #{@assembly}\n"
    out << "Total number of values: #{self.num_values}\n"
    out << "Total: #{self.total}\n"
    out << "Mean: #{self.mean}\n"
    out << "StDev: #{self.stdev}\n"
    
    # Summary for each Chromosome
  	out << self.map { |chr_id,chr_data| "Chromosome #{chr_id}: #{chr_data.length} values" }.join("\n")
      
    return out
  end
  
  # Output the data in Wig format
  def to_wig(filename)
    File.open(filename,'w') do |f|
      # TODO: should be rewritten to intelligently use step size
      f.puts Wig.track_header(@name, @description) 
          
      self.each do |chr_id,values|
        f.puts Wig.fixed_step(chr_id, values)
        f.puts values
      end
    end
  end

  # Make a Wig track header with the given name and description
  def self.track_header(name = '', description = '')
  	"track type=wiggle_0 name=\"#{name}\" description=\"#{description}\" autoScale=\"off\" visibility=\"full\""
  end
  
  # Make a Wig fixedStep chromosome header for the given ID and Chromosome
  def self.fixed_step(chr_id, chr = nil)
  	str = "fixedStep chrom=#{chr_id}"
  	
  	# fixedStep parameters
  	unless chr.nil?
  	  str << " start=#{chr.start}" if chr.start
  	  str << " step=#{chr.step}" if chr.step
	    str << " span=#{chr.span}" if chr.span
  	end
  	
  	return str
  end
  
  # Return the total number of values in this WigFile
  def num_values
  	self.collect { |chr_id,data| data.length }.sum.to_i
  end
	
	# Return the sum of all values
  def total
    self.collect { |chr_id,values| values.sum }.sum
  end

  # Return the mean of all values
  def mean
  	total.to_f / num_values
  end
	
	# Return the standard deviation of all values
  def stdev(mean = self.mean, num_values = self.num_values)
  	tss = self.map { |chr_id,values| values.to_gslv.tss(mean) }.sum
  	Math.sqrt(tss / num_values)
  end
end

class WigError < StandardError
end

##
# Base class for WigFile and BigWigFile
##
class AbstractWigFile
  include Enumerable
  
  attr_reader :name, :description, :data_file
  @@pm = Parallel::ForkManager.new(2)
  @@default_chunk_size = 200_000
  
  # Set the total number of computation processes for all Wig files (default = 2)
  def self.max_threads=(n)
    @@pm = Parallel::ForkManager.new(n.to_i)
  end
  
  # Set the default chunk size
  def self.chunk_size(n)
    @@default_chunk_size = n
  end
  
  # Open a Wig file and parse its track/chromosome information
  def initialize(filename)
    @data_file = File.expand_path(filename)
  end
  
  # Return an array of all chromosomes in this WigFile file
  def chromosomes
    raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
  end
	
	# Does this Wig file include data for chromosome chr?
	def include?(chr_id)
		raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
	end
  
  # Enumerate over the chromosomes in this Wig file
  # @DEPRECATED: Loads entire chromosomes of data, unsuitable for large genomes
  def each
    self.chromosomes.each do |chr_id| 
			yield [chr_id, chr(chr_id)]
		end
  end
  
  # Allow indexing Wig files like GenomicDatas and returning entire chromosomes
  # @DEPRECATED: Loads entire chromosomes of data, unsuitable for large genomes
  def [](chr_id)
    chr(chr_id)
  end
    
  # Load data from disk and return a Vector of values for a given chromosome
  # @DEPRECATED: Loads entire chromosomes of data, unsuitable for large genomes
  def chr(chr_id)
    raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
  end
  
  # Get the length of a chromosome from the index
  def chr_length(chr_id)
    raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
  end
  
  # Return single-bp data from the specified region
  def query(chr, start, stop)
    raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
  end
  
  # Number of values in the Wig file
  def num_values
    self.chromosomes.map { |chr| chr_length(chr) }.sum
  end
  
  # The sum of all values
  def total
    chr_totals = self.chunk_map(0) do |sum,chr,start,stop|
      sum + query(chr,start,stop).sum
    end
    
    chr_totals.sum
  end
  
  # The mean of all values
  def mean
    raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
  end
  
  # The standard deviation of all values
  def stdev(avg = self.mean)
    raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
  end
  
  # Write this Wig file to a BedGraph file
  def to_bedGraph(filename, chunk_size = 200_000)
    File.open(File.expand_path(filename), 'w') do |f|
      self.chromosomes.each do |chr|
        chunk_start = 1
        num_bases = chr_length(chr)
        while chunk_start < num_bases
          chunk_stop = [chunk_start+chunk_size-1, num_bases].min
          query(chr, chunk_start, chunk_stop).each_with_index do |value,i|
            f.puts "#{chr}\t#{chunk_start+i}\t#{chunk_start+i}\t#{value}"
          end
          chunk_start = chunk_stop + 1
        end
      end
    end
  end
  
  # Run a given block for each chromosome
  # (parallel each)
  def p_each    
    self.chromosomes.each do |chr|
      @@pm.start(chr)
      yield(chr)
      @@pm.finish(0)
    end
  end

  # Compute a given block for each chromosome
  # and return the results (parallel map)
  def p_map
    results = Array.new(self.chromosomes.length)
    
    self.chromosomes.each_with_index do |chr,i|
      @@pm.start(chr)
      results[i] = yield(chr)
      @@pm.finish(0)
    end
    
    return results
  end
  
  # Iterate over all chromosomes in chunks
  def chunk_each(chunk_size = @@default_chunk_size)
    self.chromosomes.each do |chr|
      # Run in parallel processes managed by ForkManager
      @pm.start(chr) and next 
      puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
      
      chunk_start = 1
      chr_length = self.chr_length(chr)
      while chunk_start < chr_length
        chunk_stop = [chunk_start+chunk_size-1, chr_length].min
        puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
        yield(chr, chunk_start, chunk_stop)
        chunk_start = chunk_stop + 1
      end

      @pm.finish(0)
    end
  end
  
  # Compute a given block for all chromosomes in chunks
  # and inject the results (parallel inject)
  def chunk_map(initial_value, chunk_size = @@default_chunk_size)
    results = Array.new(self.chromosomes.length)
    
    self.chromosomes.each_with_index do |chr,i|
      # Run in parallel processes managed by ForkManager
      @pm.start(chr) and next   
      puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
      
      result = initial_value
      chunk_start = 1
      chr_length = self.chr_length(chr)
      while chunk_start < chr_length
        chunk_stop = [chunk_start+chunk_size-1, chr_length].min
        puts "Computing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
        
        result = yield(result, chr, chunk_start, chunk_stop)
        
        chunk_start = chunk_stop + 1
      end
      
      results[i] = result

      @pm.finish(0)
    end
    
    return results
  end
end

##
# For documentation, see: http://genome.ucsc.edu/goldenPath/help/bigWig.html
# Analogous to WigFile, but for compressed BigWigs
##
class BigWigFile < AbstractWigFile
  attr_reader :mean, :stdev, :min, :max

  def initialize(filename)
    super(filename)
    
    output = %x[ bigWigInfo -chroms #{@data_file} ]
    raise WigError, "You must first convert your Wig file to BigWig" if output.chomp.end_with?('is not a big wig file')
    
    info = output.split("\n")
    @chromosomes = Hash.new
    info[7..-6].each do |chr_info|
      entry = chr_info.chomp.split(' ')
      chr = entry.first
      chr_size = entry.last.to_i
      @chromosomes[chr] = chr_size
    end
    
    @mean = info[-4].chomp.split(':').last.to_f
    @min = info[-3].chomp.split(':').last.to_f
    @max = info[-2].chomp.split(':').last.to_f
    @stdev = info[-1].chomp.split(':').last.to_f
  end

  # Return an array of all chromosomes in this WigFile file
  def chromosomes
    @chromosomes.keys
  end
	
	# Does this Wig file include data for chromosome chr?
	def include?(chr_id)
		@chromosomes.include?(chr_id)
	end
  
  # Load data from disk and return a Vector of values for a given chromosome
  # @DEPRECATED: Loads entire chromosomes of data, unsuitable for large genomes
  def chr(chr_id)
    chrom = Chromosome.new(chr_length(chr_id), 1, 1, 1)
    chrom[0..-1] = query(chr_id, 1, chr_length(chr_id))
    return chrom
  end
  
  # Get the length of a chromosome from the index
  def chr_length(chr_id)
    @chromosomes[chr_id]
  end
  
  # Return single-bp data from the specified region
  def query(chr, start, stop)
    %x[ bigWigSummary #{@data_file} #{chr} #{start} #{stop} #{stop-start+1} ].split(' ').map { |v| v.to_f }
  end
  
  # Write this BigWigFile to a Wig file
  def to_wig(output_file)
    puts "Converting BigWig file (#{File.basename(@data_file)}) to Wig (#{File.basename(output_file)})" if ENV['DEBUG']
    %x[ bigWigToWig #{@data_file} #{File.expand_path(output_file)} ]
  end
  
  # Output a summary about this BigWigFile
  def to_s
    str = "BigWigFile: connected to file #{@data_file}\n"
    @chromosomes.each do |chr,chr_size|
      str += "\tChromosome #{chr} (bases covered: #{chr_size})\n"
    end
    
    str += "Mean:\t#{@mean}"
    str += "Standard deviation:\t#{@stdev}"
    str += "Min:\t#{@min}"
    str += "Max:\t#{@max}"
    
    return str
  end
end

##
# Lazy-load a Wig file, reading data into memory by chromosome
# only as needed to conserve memory
##
class WigFile < AbstractWigFile
  # Open a Wig file and parse its track/chromosome information
  def initialize(filename)
    super(filename)
    
    # Load the track information from the first line
    File.open(@data_file) do |f|
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
		
		# Find chromosome header lines and index (ghetto B-tree index)
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
  
  # Return an array of all chromosomes in this WigFile file
  def chromosomes
    @index.keys
  end
	
	# Does this Wig file include data for chromosome chr?
	def include?(chr_id)
		@index.include?(chr_id)
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
  
  # The mean of all values
  def mean
    total / num_values
  end
  
  # The standard deviation of all values
  def stdev(avg = self.mean)
    chr_deviances = self.chunk_map(0) do |sum,chr,start,stop|
      sum + query(chr,start,stop).map { |elem| (elem-avg)**2 }.sum
    end
    
    Math.sqrt(chr_deviances.sum / num_values)
  end
  
	# Output a summary about this WigFile
  def to_s
    str = "WigFile: connected to file #{@data_file}\n"
    @index.each do |chr,line|
      str += "\tChromosome #{chr} (lines: #{line}..#{chr_stop(chr)})\n"
    end
    
    return str
  end
  
  # Convert this WigFile to a BigWigFile
  def to_bigwig(output_file, assembly)
    puts "Converting Wig file (#{File.basename(@data_file)}) to BigWig (#{File.basename(output_file)})" if ENV['DEBUG']
    %x[ wigToBigWig #{@data_file} #{assembly.len_file} #{output_file} ]
  end
end
