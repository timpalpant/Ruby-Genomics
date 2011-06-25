require 'enumerator'
require 'unix_file_utils'
require 'parallelizer'
require 'contig'
require 'assembly'
require 'stringio'
require 'wig_math'
require 'ucsc_tools'

class WigError < StandardError
end

##
# Base class for WigFile and BigWigFile
##
class AbstractWigFile
  include Enumerable
  include WigMath
  
  attr_reader :track_header, :data_file
  
  ##
  # INSTANCE METHODS
  ##
  
  # Open a Wig file and parse its track/chromosome information
  def initialize(filename)
    @data_file = File.expand_path(filename)
    @track_header = UCSCTrackHeader.new(:type => 'wiggle_0')
    @contigs_index = Array.new
  end
  
  # Open a Wig file with an optional block
  def open(filename, &block)
    wig = self.new(filename)
  
    if block
      yield wig
    else
      return wig
    end
    
    # Place any cleanup operations here
  end
  
  # Enumerate over the contigs in this Wig file
  def each
    @contigs_index.each do |contig_info|
      yield query(contig_info.chr, contig_info.start, contig_info.stop)
    end
  end
  
  # Enumerate over chunks in this Wig file
  # NOTE: The order of the chunks is not guaranteed
  def each_chunk
    @contigs_index.each do |contig_info|
      chunk_start = contig_info.start
      while chunk_start <= contig_info.stop
        chunk_stop = [chunk_start+200_000-1, contig_info.stop].min
        puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
        
        yield query(contig_info.chr, chunk_start, chunk_stop)

        chunk_start = chunk_stop + 1
      end
    end
  end
  
  # Return an array of all chromosomes in this Wig file
  def chromosomes
    @contig_index.collect { |contig_info| contig_info.chr }.uniq
  end
  
  # Does this Wig file include data for a given locus?
  def include?(chr, start = nil, stop = nil)
    if start.nil?
      chromosomes.include?(chr)
    elsif stop.nil?
      chromsomes.include?(chr) and start >= chr_start(chr)
    else
      chromsomes.include?(chr) and start >= chr_start(chr) and stop <= chr_stop(chr)
    end
  end
  
  # Get the contigs for a chromosome
  def chr_contigs(query_chr)
    raise WigError, "Wig does not include data for chromosome #{query_chr}" unless include?(query_chr)
    return @contigs_index.select { |contig_info| contig_info.chr == query_chr }
  end
  
  # Get the lowest base pair of data for a chromosome
  def chr_start(query_chr)
    chr_contigs(query_chr).collect { |contig_info| contig_info.start }.min
  end
  
  # Get the highest base pair of data for a chromosome
  def chr_stop(chr)
    chr_contigs(query_chr).collect { |contig_info| contig_info.stop }.max
  end
  
  # Compute some transformation on the values in this file, and output the result to disk
  def transform(output_file)
    # Write the output file header
    header_file = output_file+'.header'
    File.open(header_file, 'w') do |f|
      f.puts @track_header
    end

    # Keep track of all the temporary intermediate files (header first)
    tmp_files = [header_file]
    @contigs_index.each { |contig_info| tmp_files << output_file+'.'+contig_info.chr+'.'+contig_info.start }
  
    begin
      # Iterate by contig
      @contigs_index.p_each do |contig_info|
        puts "\nProcessing contig #{contig_info}" if ENV['DEBUG']

        # Write the header
        chr_temp_file = output_file+'.'+contig_info.chr+'.'+contig_info.start
        File.open(chr_temp_file, 'w') do |f|
          #f.puts Contig.new(0, chr, 1, 1, 1).to_s
        end
        
        chunk_start = contig_info.start
        while chunk_start <= contig_info.stop
          chunk_stop = [chunk_start+200_000-1, contig_info.stop].min
          puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
          
          output = yield(contig_info.chr, chunk_start, chunk_stop)
          
          # Write this chunk to disk
          File.open(chr_temp_file, 'a') do |f|
            f.puts output.map { |value| value.to_s(5) }.join("\n")
          end
          
          chunk_start = chunk_stop + 1
        end
      end
      
      # Concatenate all of the temp file pieces into the final output
      File.cat(tmp_files, output_file)
    rescue
      raise WigError, "Error transforming Wig file!"
    ensure
      # Delete the individual temp files created by each process
      tmp_files.each { |filename| File.delete(filename) if File.exist?(filename) }
    end
  end
  
  # Output a summary about this BigWigFile
  def to_s
    str = StringIO.new
    str << "AbstractWigFile: connected to file #{@data_file}\n"
    @contigs_index.each do |chr,start,stop|
      str << "\tContig - #{chr}:#{start}-#{stop}\n"
    end
    
    str << "Mean:\t#{mean}\n"
    str << "Standard deviation:\t#{stdev}"
    
    return str.string
  end
  
  ##
  # ABSTRACT METHODS
  ##
  
  # Return a Contig of data from the specified region
  def query(chr, start, stop)
    raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
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
    
    contig_info = ContigInfo.new
    info[7..-6].each do |line|
      entry = line.chomp.split(' ')
      contig_info.chr = entry.first
      # Assume that all chromosomes start at base 1 until
      # we find a way to get the first base
      contig_info.start = 1
      contig_info.stop = entry.last.to_i
      @contigs_index << contig_info
    end
    
    # bigWigInfo doesn't calculate mean/stdev accurately enough when values are small (10^-7)
    #@mean = info[-4].chomp.split(':').last.to_f
    #@stdev = info[-1].chomp.split(':').last.to_f
    @min = info[-3].chomp.split(':').last.to_f
    @max = info[-2].chomp.split(':').last.to_f
  end
  
  # Convert the output of #transform back to BigWig
  alias :super_transform :transform
  def transform(output_file, assembly, &block)
    super_transform(output_file, &block)
    
    # Convert the output Wig file to BigWig
    tmp_file = output_file + '.tmp'
    begin
      WigFile.to_bigwig(output_file, tmp_file, assembly)
    
      # Delete the temporary intermediate Wig file by moving the BigWig on top of it
      FileUtils.move(tmp_file, output_file)
    rescue
      # Cleanup the files if something went wrong
      File.delete(tmp_file) if File.exist?(tmp_file)
      File.delete(output_file) if File.exist?(output_file)
    end
  end
  
  ##
  # QUERY METHODS
  ##
  
  # Return a Contig of data from the specified region
  def query(chr, start = nil, stop = nil)
    # Don't query off the ends of chromosomes
    start = chr_start(chr) if start.nil?
    stop = chr_stop(chr) if stop.nil?
    raise WigError, "BigWig does not contain data for the interval #{chr}:#{start}-#{stop}" if not include?(chr, start, stop)

    # Data is 0-indexed
    # bigWigSummary segfaults if query is too big, so use nice, bite-sized chunks
    query_start = start-1
    contig = Contig.new(chr)
    while query_start <= stop-1
      query_stop = [query_start+200_000-1, stop-1].min
      num_values = query_stop-query_start+1
      chunk = %x[ bigWigSummary #{@data_file} #{chr} #{query_start} #{query_stop} #{num_values} 2>&1 ]
      raise WigError, "BigWig does not contain data for the interval #{chr}:#{query_start+1}-#{query_stop+1}" if chunk.start_with?('no data in region')
      
      # Store the chunk of values in the Contig
      count = 0
      chunk.split(' ').each_with_index do |value,i|
        contig.set(query_start+i, value.to_f)
        count += 1
      end
      raise WigError, "BigWig query did not return the expected number of values! (#{count} != #{num_values})" if count != num_values
      
      query_start = query_stop + 1
    end
    
    return contig
  end

  # Return the average value for the specified region
  def query_average(chr, start, stop)
    # Don't query off the ends of chromosomes
    raise WigError, "BigWig does not contain data for the interval #{chr}:#{start}-#{stop}" if not include?(chr, start, stop)

    %x[ bigWigSummary #{@data_file} #{chr} #{start-1} #{stop-1} 1 ].to_f
  end
  
  ##
  # OUTPUT METHODS
  ##
  
  # Write a BigWigFile to a Wig file
  def self.to_wig(input_file, output_file)
    puts "Converting BigWig file (#{File.basename(input_file)}) to Wig (#{File.basename(output_file)})" if ENV['DEBUG']
    
    header_file = File.expand_path(output_file + '.header')
    data_file = File.expand_path(output_file + '.data')
    
    begin
      # Write a track header
      File.open(header_file, 'w') do |f|
        f.puts UCSCTrackHeader.new(:name => File.basename(output_file))
      end
      
      # Extract the data with UCSC tools
      %x[ bigWigToWig #{input_file} #{data_file} ]
      
      # Cat the two parts together
      File.cat([header_file, data_file], output_file)
    rescue
      raise WigError, "Error converting BigWig file to Wig"
    ensure
      # Delete the two temp files
      File.delete(header_file) if File.exist?(header_file)
      File.delete(data_file) if File.exist?(data_file)
    end
  end

  # Write this BigWig to a BedGraph
  def self.to_bedgraph(input_file, output_file)
    puts "Converting BigWig file (#{File.basename(input_file)}) to BedGraph (#{File.basename(output_file)})" if ENV['DEBUG']
    %x[ bigWigToBedGraph #{input_file} #{output_file} ]
  end
end


##
# An ASCII text Wiggle file
##
class WigFile < AbstractWigFile
  # Open a Wig file and parse its track/contig information
  def initialize(filename)
    super(filename)
    
    # Load the track information from the first line
    File.open(@data_file) do |f|
      begin
        @track_header = UCSCTrackHeader.parse(f.gets)
      rescue UCSCTrackHeaderError
        puts "Error parsing track header" if ENV['DEBUG']
      end
    end
    
    # Call grep to load the contig information
    puts 'Indexing Contig header lines' if ENV['DEBUG']
    File.grep_with_linenum(@data_file, 'chrom') do |line_num,line|      
      begin
        info = ContigInfo.parse(line)
        info.line_start = line_num
        @contigs_index << info
      rescue
        puts "Not a valid fixedStep/variableStep header" if ENV['DEBUG']
        next
      end
    end
    
    # Now find the start and stop of each Contig
    @contigs_index.each do |contig_info|
      # fixedStep lines give the start, so we just need to find the stop
      if contig_info.type == ContigInfo::FIXED_STEP
        # Get the line number of the next contig in the file
        next_contig_line = @contigs_index.select { |info| info.line_start > contig_info.line_start }.collect { |info| info.line_start }.sort.first
        offset = 0
        last_line = String.new
        while last_line.empty?
          offset += 1
          last_line = File.lines(@data_file, next_contig_line-offset, next_contig_line-offset).first.chomp
        end
        contig_info.line_stop = next_contig_line - offset
        
        # Calculate the stop based on the number of values and the step/span sizes
        num_values = contig_info.line_stop - contig_info.line_start
        contig_info.stop = contig_info.step*(num_values-1) + contig_info.span
      # for variableStep lines, we need to find the start and stop
      else
        # Find the start, i.e. the first base pair with data
        first_line = String.new
        offset = 0
        while first_line.empty?
          offset += 1
          first_line = File.lines(@data_file, contig_info.line_start+offset, contig_info.line_start+offset).first
        end
        contig_info.start = first_line.split("\t").first.to_i
        
        # Find the stop, i.e. the last base pair with data
        # Get the line number of the next contig in the file
        next_contig_line = @contigs_index.select { |info| info.line_start > contig_info.line_start }.collect { |info| info.line_start }.sort.first
        # Get the line immediately before it (with data)
        offset = 0
        last_line = String.new
        while last_line.empty?
          offset += 1
          last_line = File.lines(@data_file, next_contig_line-offset, next_contig_line-offset).first.chomp
        end
        contig_info.line_stop = next_contig_line - offset
        contig_info.stop = last_line.split("\t").first.to_i + contig_info.span
      end
    end
  
    # Raise an error if no chromosomes were found
    raise WigError, "No fixedStep/variableStep headers found in Wig file!" if @contigs_index.length == 0
  end
  
  ##
  # QUERY METHODS
  ##
  
  # Return single-bp data from the specified region
  def query(chr, start = nil, stop = nil)
    start = chr_start(chr) if start.nil?
    stop = chr_stop(chr) if stop.nil?
    raise WigError, "Chromosome #{chr} not found in Wig file #{@data_file}!" unless include?(chr, start, stop)
    
    output = Contig.new(chr)
    
    # Find the relevant contigs for the requested interval
    @contigs_index.select { |info| info.chr == chr and (info.stop >= start and info.start <= stop) }.sort_by { |info| info.start }.each do |info|
      # Figure out which bases are covered by this Contig
      low = [start, info.start].max
      high = [stop, info.stop].min
      
      if info.type == ContigInfo::FIXED_STEP
        # Figure out what lines in the file we need to get those bases
        # (info.start_line + 1) is the first line of data
        start_line = (info.start_line + 1) + (low-info.start)/info.step
        stop_line = start_line + (high-low)/info.step
        
        # Query the file for the lines and store them in the Contig
        bp = info.start + (low-info.start)/info.step
        File.lines(@data_file, start_line, stop_line) do |line|
          value = line.to_f
          start_base = [bp, low].max
          stop_base = [bp+info.span-1, stop].min
          (start_base..stop_base).each { |base| output.set(bp, value) }
          bp += info.step
        end
      else
        # With variableStep, there's no way to know exactly what lines are needed
        # so iterate through them all until we find what we want
        File.lines(@data_file, info.line_start, info.line_stop) do |line|
          entry = line.split("\t")
          start_base = entry.first.to_i
          value = entry.last.to_f
          stop_base = [bp+info.span-1, high].min
          
          if start_base >= low
            (start_base..stop_base).each { |bp| output.set(bp, value) }
          end
          
          # Don't need to go further in the Contig if we've reached our last base
          break if stop_base == high
        end
      end
    end
    
    return output
  end
  
  ##
  # OUTPUT METHODS
  ##
  
  # Convert this WigFile to a BigWigFile
  def to_bigwig(output_file, assembly)
    WigFile.to_bigwig(@datafile, output_file, assembly)
  end

  # Convert this WigFile to a BedGraph
  def to_bedgraph(output_file)
    WigFile.to_bedgraph(@datafile, output_file)
  end
  
  # For converting wigs to BigWigs without having to load (index them) first
  def self.to_bigwig(input_file, output_file, assembly)
    puts "Converting Wig file (#{File.basename(input_file)}) to BigWig (#{File.basename(output_file)})" if ENV['DEBUG']
    %x[ wigToBigWig -clip #{File.expand_path(input_file)} #{File.expand_path(assembly.len_file)} #{File.expand_path(output_file)} ]
  end
  
  # For converting wigs to BedGraph without having to load (index them) first
  # Also creates the most compact BedGraph possible by joining equal neighbors
  def self.to_bedgraph(input_file, output_file)
    puts "Converting Wig file #{File.basename(input_file)} to BedGraph #{File.basename(output_file)}" if ENV['DEBUG']
    raise "Not yet implemented"
  end
end

##
# Holds info about a Contig in a WigFile
##
class ContigInfo
  attr_accessor :type, :chr, :start, :stop, :step, :span, :line_start, :line_stop
  
  FIXED_STEP = 'fixedStep'
  VARIABLE_STEP = 'variableStep'
  
  def initialize(type = FIXED_STEP, chr = 'unknown', start = 1, stop = 1, step = 1, span = 1, line_start = nil, line_stop = nil)
    @type = type
    @chr = chr
    @start = start
    @stop = stop
    @step = step
    @span = span
    @line_start= line_start
    @line_stop = line_stop
  end
  
  # Parse a fixedStep/variableStep line
  def self.parse(line)
    # Remove any leading/trailing whitespace
    line.chomp!
    
    # Store the type of Contig (fixedStep / variableStep)
    info = self.new
    if line.start_with?(FIXED_STEP)
      info.type = FIXED_STEP
    elsif line.start_with?(VARIABLE_STEP)
      info.type = VARIABLE_STEP
    else
      raise WigError, "Not a valid fixedStep/variableStep line!"
    end
  
    # Parse the other tokens
    line.split(' ').each do |opt|
      keypair = opt.split('=')
      key = keypair.first
      value = keypair.last
      
      case key
        when 'chrom'
          info.chr = value
        when 'start'
          info.start = value.to_i
        when 'step'
          info.step = value.to_i
        when 'span'
          info.span = value.to_i
      end
    end

    return info
  end
end

class WigError < StandardError
end
