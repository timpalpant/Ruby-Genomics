require 'enumerator'
require 'tmpdir'
require 'unix_file_utils'
require 'contig'
require 'genomic_index_error'
require 'stats'
require 'forkmanager'
require 'assembly'
require 'stringio'

##
# Shouldn't ever load a Wig file completely into memory, so just methods to convert and create Wig files
##
class Wig
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

  # For converting wigs to BigWigs without having to load (index them) first
  def self.to_bigwig(input_file, output_file, assembly)
    puts "Converting Wig file (#{File.basename(input_file)}) to BigWig (#{File.basename(output_file)})" if ENV['DEBUG']
    %x[ wigToBigWig -clip #{File.expand_path(input_file)} #{File.expand_path(assembly.len_file)} #{File.expand_path(output_file)} ]
  end
  
  # For converting wigs to BedGraph without having to load (index them) first
  # Also creates the most compact BedGraph possible by joining equal neighbors
  def self.to_bedGraph(input_file, output_file)
    puts "Converting Wig file #{File.basename(input_file)} to BedGraph #{File.basename(output_file)}" if ENV['DEBUG']
    File.open(File.expand_path(output_file), 'w') do |f|
      chr, current, step, span = nil, 1, 1, 1
      start, prev_value = 1, 0
      File.foreach(File.expand_path(input_file)) do |line|
        if line.start_with?('track')
          next
        elsif line.start_with?('variable')
          raise WigError, "Only fixedStep-format Wig files are supported at this time"
        elsif line.start_with?('fixed')
          line.split(' ').each do |opt|
            keypair = opt.split('=')
            key = keypair.first
            value = keypair.last
            
            if key == 'chrom'
              chr = value
              puts "Processing chromosome #{chr}" if ENV['DEBUG']
            elsif key == 'start'
              current = value.to_i
            elsif key == 'step'
              step = value.to_i
            elsif key == 'span'
              span = value.to_i
            end
          end        
        else
          value = line.chomp.to_f
          if value == prev_value
            current += step
          else
            f.puts "#{chr}\t#{start}\t#{current+span-1}\t#{prev_value}"
            current += step
            start = current
            prev_value = value
          end
        end
      end
    end
  end
end

##
# Shouldn't load BigWig files into memory since they can be randomly accessed
# so just methods to convert them into other things
##
class BigWig
  # Write a BigWigFile to a Wig file
  def self.to_wig(input_file, output_file)
    puts "Converting BigWig file (#{File.basename(input_file)}) to Wig (#{File.basename(output_file)})" if ENV['DEBUG']
    %x[ bigWigToWig #{input_file} #{File.expand_path(output_file)} ]
  end

  # Write this BigWig to a BedGraph
  def self.to_bedGraph(input_file, output_file)
  puts "Converting BigWig file (#{File.basename(input_file)}) to BedGraph (#{File.basename(output_file)})" if ENV['DEBUG']
  %x[ bigWigToBedGraph #{input_file} #{output_file} ]
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
  @@pm = Parallel::ForkManager.new(2, {'tempdir' => Dir.tmpdir})
  @@default_chunk_size = 200_000
  
  # Set the total number of computation processes for all Wig files (default = 2)
  def self.max_threads=(n)
    @@pm = Parallel::ForkManager.new(n.to_i, {'tempdir' => Dir.tmpdir})
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

    return chr_totals.sum
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
  
  # Run a given block for each chromosome
  # (parallel each)
  def p_each    
    self.chromosomes.each do |chr|
      @@pm.start(chr)
      yield(chr)
      @@pm.finish(0)
    end

    @@pm.wait_all_children
  end

  # Compute a given block for each chromosome
  # and return the results (parallel map)
  def p_map
    results = Array.new(self.chromosomes.length)
    
    @@pm.run_on_finish do |pid,exit_code,id,exit_signal,core_dump,data|
      results[id] = data['output']
    end

    self.chromosomes.each_with_index do |chr,i|
      @@pm.start(chr)
      result = yield(chr)
      @@pm.finish(0, {'output' => result})
    end
    
    @@pm.wait_all_children

    return results
  end
  
  # Iterate over all chromosomes in chunks
  def chunk_each(chunk_size = @@default_chunk_size)
    self.chromosomes.each do |chr|
      # Run in parallel processes managed by ForkManager
      @@pm.start(chr) and next 
      puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
      
      chunk_start = 1
      chr_length = self.chr_length(chr)
      while chunk_start < chr_length
        chunk_stop = [chunk_start+chunk_size-1, chr_length].min
        puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
        yield(chr, chunk_start, chunk_stop)
        chunk_start = chunk_stop + 1
      end

      @@pm.finish(0)
    end

    @@pm.wait_all_children
  end
  
  # Compute a given block for all chromosomes in chunks
  # and inject the results (parallel inject)
  def chunk_map(initial_value, chunk_size = @@default_chunk_size)
    results = Array.new(self.chromosomes.length)
    
    @@pm.run_on_finish do |pid,exit_code,id,exit_signal,core_dump,data|
      results[id] = data['output']
    end

    self.chromosomes.each_with_index do |chr,i|
      # Run in parallel processes managed by ForkManager
      @@pm.start(i) and next   
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
      
      @@pm.finish(0, {'output' => result})
    end

    @@pm.wait_all_children
    
    return results
  end
end

##
# For documentation, see: http://genome.ucsc.edu/goldenPath/help/bigWig.html
# Analogous to WigFile, but for compressed BigWigs
##
class BigWigFile < AbstractWigFile
  attr_reader :min, :max

  def initialize(filename)
    super(filename)
    
    info = %x[ bigWigInfo -chroms #{@data_file} ].split("\n")
    raise WigError, "You must first convert your Wig file to BigWig" if info.length < 8
    
    @chromosomes = Hash.new
    info[7..-6].each do |chr_info|
      entry = chr_info.chomp.split(' ')
      chr = entry.first
      chr_size = entry.last.to_i
      @chromosomes[chr] = chr_size
    end
    
    # bigWigInfo doesn't calculate mean/stdev accurately enough when values are small (10^-7)
    #@mean = info[-4].chomp.split(':').last.to_f
    #@stdev = info[-1].chomp.split(':').last.to_f
    @min = info[-3].chomp.split(':').last.to_f
    @max = info[-2].chomp.split(':').last.to_f
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
    query(chr_id, 1, chr_length(chr_id))
  end
  
  # Get the length of a chromosome from the index
  def chr_length(chr_id)
    @chromosomes[chr_id]
  end
  
  # Return single-bp data from the specified region
  def query(chr, start, stop)
    # Don't query off the ends of chromosomes
    raise GenomicIndexError, "BigWig does not contain data for the interval #{chr}:#{start}-#{stop}" if not include?(chr) or start < 1 or stop > chr_length(chr)

    # Allow Crick queries
    low = [start, stop].min
    high = [start, stop].max

    # Data is 0-indexed
    # bigWigSummary segfaults if query is too big, so use nice bite-size chunks
    query_start = low-1
    result = StringIO.new
    while query_start <= high-1
      query_stop = [query_start+@@default_chunk_size-1, high-1].min
      chunk = %x[ bigWigSummary #{@data_file} #{chr} #{query_start} #{query_stop} #{query_stop-query_start+1} 2>&1 ]
      raise GenomicIndexError, "BigWig does not contain data for the interval #{chr}:#{query_start+1}-#{query_stop+1}" if chunk.start_with?('no data in region')
      result << ' ' << chunk
      query_start = query_stop + 1
    end

    values = result.string.split(' ').map { |v| v.to_f }
    values.reverse! if start > stop
    return values.to_contig(start, 1, 1)
  end

  # Return the average value for the specified region
  def query_average(chr, start, stop)
    # Don't query off the ends of chromosomes
    raise GenomicIndexError, "BigWig does not contain data for the interval #{chr}:#{start}-#{stop}" if start < 1 or stop > chr_length(chr)

    low = [start, stop].min
    high = [start, stop].max

    %x[ bigWigSummary #{@data_file} #{chr} #{low-1} #{high-1} 1 ].to_f
  end
  
  # Output a summary about this BigWigFile
  def to_s
    str = "BigWigFile: connected to file #{@data_file}\n"
    @chromosomes.each do |chr,chr_size|
      str += "\t#{chr} (bases covered: #{chr_size})\n"
    end
    
    str += "Mean:\t#{mean}\n"
    str += "Standard deviation:\t#{stdev}\n"
    str += "Min:\t#{@min}\n"
    str += "Max:\t#{@max}\n"
    
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
      #raise WigError, "Only fixedStep-format Wig files are supported at this time" if header_line.start_with?('variable')
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
    raise WigError, "No chromosome fixedStep headers found in Wig file!" if @index.length == 0
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
    return Contig.load_wig(@data_file, chr_start(chr_id), chr_stop(chr_id))
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
    parsed = Contig.parse_wig_header(header)
    
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
    parsed[0..-1] = values
    
    return parsed
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
    Wig.to_bigwig(@datafile, output_file, assembly)
  end

  # Convert this WigFile to a BedGraph
  def to_bedgraph(output_file)
    Wig.to_bedGraph(@datafile, output_file)
  end
end
