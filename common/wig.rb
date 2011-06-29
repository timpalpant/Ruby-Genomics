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
  
  # Put any cleanup operations here
  def close
  end
  
  # Open a Wig file with an optional block
  def self.open(filename, &block)
    wig = self.new(filename)
  
    if block
      yield wig
    else
      return wig
    end
    
    wig.close
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
    @contigs_index.collect { |contig_info| contig_info.chr }.uniq
  end
  
  # Does this Wig file include data for a given locus?
  def include?(chr, start = nil, stop = nil)
    if start.nil?
      chromosomes.include?(chr)
    elsif stop.nil?
      chromosomes.include?(chr) and start >= chr_start(chr) and start <= chr_stop(chr)
    else
      chromosomes.include?(chr) and start >= chr_start(chr) and stop <= chr_stop(chr)
    end
  end
  
  # Get the lowest base pair of data for a chromosome
  def chr_start(query_chr)
    chr_contigs(query_chr).collect { |contig_info| contig_info.start }.min
  end
  
  # Get the highest base pair of data for a chromosome
  def chr_stop(query_chr)
    chr_contigs(query_chr).collect { |contig_info| contig_info.stop }.max
  end
  
  # Compute some transformation on the values in this file, and output the result to disk
  def transform(output_file)
    # Write the output file header
    header_file = output_file+'.header'
    File.open(header_file, 'w') do |f|
      f.puts @track_header.to_s
    end

    # Keep track of all the temporary intermediate files (header first)
    tmp_files = [header_file]
    @contigs_index.each { |contig_info| tmp_files << "#{output_file}.#{contig_info.chr}.#{contig_info.start}" }
  
    begin
      # Iterate by contig
      @contigs_index.p_each do |contig_info|
        puts "\nProcessing contig #{contig_info}" if ENV['DEBUG']

        # Write the header
        chr_temp_file = "#{output_file}.#{contig_info.chr}.#{contig_info.start}"
        File.open(chr_temp_file, 'w') do |f|
          f.puts "fixedStep chrom=#{contig_info.chr} start=#{contig_info.start} step=1 span=1"
        end
        
        chunk_start = contig_info.start
        while chunk_start <= contig_info.stop
          chunk_stop = [chunk_start+200_000-1, contig_info.stop].min
          puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
          
          output = yield(contig_info.chr, chunk_start, chunk_stop)
          
          # Write this chunk to disk
          File.open(chr_temp_file, 'a') do |f|
            f.puts output.map { |value| value ? value.to_s(5) : 'NaN' }.join("\n")
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
    str << "Wig: connected to file #{File.basename(@data_file)}\n"
    @contigs_index.each { |contig_info| str << "#{contig_info}\n" }
    str << "Mean:\t#{mean}\n"
    str << "Standard deviation:\t#{stdev}"
    
    return str.string
  end
  
  def inspect
    to_s
  end
  
  ##
  # ABSTRACT METHODS
  ##
  
  # Return a Contig of data from the specified region
  def query(chr, start, stop)
    raise WigError, "Should be overridden in a base class (BigWigFile/WigFile)!"
  end
  
  ##
  # HELPER METHODS
  ##
  private
  
  # Get the contigs for a chromosome
  def chr_contigs(query_chr)
    raise WigError, "Wig does not include data for chromosome #{query_chr}" unless include?(query_chr)
    return @contigs_index.select { |contig_info| contig_info.chr == query_chr }
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
    
    info[7..-6].each do |line|
      entry = line.chomp.split(' ')
      
      contig_info = ContigInfo.new
      contig_info.chr = entry.first
      chr_length = entry.last.to_i
      contig_info.start = find_start_base(contig_info.chr, chr_length)
      contig_info.stop = find_stop_base(contig_info.chr, chr_length)
      
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
  def query(chr, start = nil, stop = nil, type = 'mean')
    # Don't query off the ends of chromosomes
    start = chr_start(chr) if start.nil?
    stop = chr_stop(chr) if stop.nil?
    raise WigError, "BigWig does not contain data for the interval #{chr}:#{start}-#{stop}" if not include?(chr, start, stop)

    # Data is 0-indexed and half-open
    # bigWigSummary segfaults if query is too big, so use nice, bite-sized chunks
    query_start = start-1
    contig = Contig.new(chr)
    while query_start <= stop
      query_stop = [query_start+200_000, stop].min
      num_values = query_stop-query_start
      chunk = %x[ bigWigSummary -type=#{type} #{@data_file} #{chr} #{query_start} #{query_stop} #{num_values} 2>&1 ]
      raise WigError, "BigWig does not contain data for the interval #{chr}:#{query_start+1}-#{query_stop+1}" if chunk.start_with?('no data in region')
      
      # Store the chunk of values in the Contig
      count = 0
      chunk.split(' ').each_with_index do |value,i|
        count += 1
        next if value == 'n/a' or value == 'nan'
        contig.set(query_start+i+1, value.to_f)
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
    %x[ bigWigSummary #{@data_file} #{chr} #{start-1} #{stop} 1 ].to_f
  end
  
  ##
  # OUTPUT METHODS
  ##
  
  def to_wig(output_file)
    BigWig.to_wig(@data_file, output_file)
  end
  
  def to_bedgraph(output_file)
    BiGWig.to_bedgraph(@data_file, output_file)
  end
  
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

  ##
  # HELPER METHODS
  ##

  private

  # Find the first base pair with data
  def find_start_base(chr, chr_length)
    # Start at the first base pair and look forwards
    chunk_size = 200_000
    bp = 1
    while bp <= chr_length
      start = bp-1
      stop = [start+chunk_size, chr_length-1].min
      num_values = stop - start
      result = %x[ bigWigSummary -type=coverage #{@data_file} #{chr} #{start} #{stop} #{num_values} 2>&1 ]
      if not result.start_with?('no data')
        start = result.split(' ').find_index('1')
        return bp+start unless start.nil?
      end
      
      bp += chunk_size
    end
    
    raise WigError, "Could not find start base pair in BigWig file #{File.basename(@data_file)} for chromosome #{chr}"
  end

  # Find the last base pair with data
  def find_stop_base(chr, chr_length)
    # Start at the last base pair and look backwards
    chunk_size = 200_000
    bp = chr_length
    while bp >= 1
      start = [1, bp-chunk_size-1].max
      stop = bp-1
      num_values = stop - start
      result = %x[ bigWigSummary -type=coverage #{@data_file} #{chr} #{start} #{stop} #{num_values} 2>&1 ]
      if not result.start_with?('no data')
        start = result.split(' ').reverse.find_index('1')
        return bp-start-1 unless start.nil?
      end
      
      bp -= chunk_size
    end
    
    raise WigError, "Could not find stop base pair in BigWig file #{File.basename(@data_file)} for chromosome #{chr}"
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
        info.header_line = line_num
        line_start = line_num
        first_line = String.new
        while first_line.nil? or first_line.chomp.empty?
          line_start += 1
          first_line = File.lines(@data_file, line_start, line_start).first
        end
        info.line_start = line_start
        @contigs_index << info
      rescue
        puts "Not a valid fixedStep/variableStep header" if ENV['DEBUG']
      end
    end
    
    # Now find the start and stop of each Contig
    @contigs_index.each do |contig_info|
      # fixedStep lines give the start, so we just need to find the stop
      if contig_info.type == ContigInfo::FIXED_STEP
        # Get the line number of the next contig in the file
        next_contig_line = @contigs_index.select { |info| info.header_line > contig_info.header_line }.collect { |info| info.header_line }.sort.first
        
        # Get the line immediately before it (with data)
        offset = 0
        
        # If there is no contig after this one, then read to EOF
        if next_contig_line.nil?
          next_contig_line = File.num_lines(@data_file)
          offset -= 1
        end
        
        last_line = String.new
        while last_line.nil? or last_line.chomp.empty?
          offset += 1
          last_line = File.lines(@data_file, next_contig_line-offset, next_contig_line-offset).first
        end
        contig_info.line_stop = next_contig_line - offset
        
        # Calculate the stop based on the number of values and the step/span sizes
        num_values = contig_info.line_stop - contig_info.line_start
        contig_info.stop = contig_info.start + contig_info.step*(num_values) + contig_info.span - 1
      # for variableStep lines, we need to find the start and stop
      else
        # Find the start, i.e. the first base pair with data
        contig_info.start = File.lines(@data_file, contig_info.line_start, contig_info.line_start).first.chomp.split("\t").first.to_i
        
        # Find the stop, i.e. the last base pair with data
        # Get the line number of the next contig in the file
        next_contig_line = @contigs_index.select { |info| info.header_line > contig_info.header_line }.collect { |info| info.header_line }.sort.first
        
        # Get the line immediately before it (with data)
        offset = 0
        
        # If there is no contig after this one, then read to EOF
        if next_contig_line.nil?
          next_contig_line = File.num_lines(@data_file)
          offset -= 1
        end
        
        last_line = String.new
        while last_line.nil? or last_line.chomp.empty?
          offset += 1
          last_line = File.lines(@data_file, next_contig_line-offset, next_contig_line-offset).first
        end
        contig_info.line_stop = next_contig_line - offset
        contig_info.stop = last_line.split("\t").first.to_i + contig_info.span - 1
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
      # Clamp to bases that are covered by this Contig
      low = [start, info.start].max
      high = [stop, info.stop].min
      
      if info.type == ContigInfo::FIXED_STEP
        # Figure out what lines in the file we need to get those bases
        start_line = info.line_for_bp(low)
        stop_line = info.line_for_bp(high)
        
        # Query the file for the lines and store them in the Contig
        bp = info.bp_for_line(start_line)
        File.lines(@data_file, start_line, stop_line) do |line|
          value = line.to_f
          (bp..bp+info.span-1).each { |base| output.set(base, value) if low <= base and base <= high }
          bp += info.step
        end
      else
        # With variableStep, there's no way to know exactly what lines are needed
        # so iterate through them all until we find what we want
        File.lines(@data_file, info.line_start, info.line_stop) do |line|
          entry = line.split("\t")
          start_base = entry.first.to_i
          value = entry.last.to_f
          stop_base = [start_base+info.span-1, high].min
          
          if start_base >= low
            (start_base..stop_base).each { |base| output.set(base, value) }
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
  attr_accessor :type, :chr, :start, :stop, :step, :span, :header_line, :line_start, :line_stop
  
  FIXED_STEP = 'fixedStep'
  VARIABLE_STEP = 'variableStep'
  
  def initialize(type = FIXED_STEP, chr = 'unknown', start = 1, stop = 1, step = 1, span = 1, header_line = nil, line_start = nil, line_stop = nil)
    @type = type
    @chr = chr
    @start = start
    @stop = stop
    @step = step
    @span = span
    @header_line = header_line
    @line_start= line_start
    @line_stop = line_stop
  end
  
  # Return which line contains data for a given base pair
  def line_for_bp(bp)
    raise WigError, "Cannot compute the line for a base pair with variableStep contigs" if @type == VARIABLE_STEP
    raise WigError, "Contig does not contain data for base pair #{bp}" if bp < start or bp > stop
    @line_start + (bp-@start)/@step
  end
  
  # Return the start base pair of a given line number
  def bp_for_line(line_num)
    raise WigError, "Cannot compute the base pair for a line with variableStep contigs" if @type == VARIABLE_STEP
    raise WigError, "Contig does not include line number #{line_num}" if line_num < @line_start or line_num > @line_stop
    @start + @step * (line_num-@line_start)
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

  def to_s
    s = StringIO.new
    s << "#{@type}, chrom=#{@chr}, start=#{@start}, stop=#{@stop}"
    s << ", header_line=#{@header_line}" if @header_line
    s << ", start_line=#{@line_start}" if @line_start
    s << ", stop_line=#{@line_stop}" if @line_stop
    return s.string
  end
end

class WigError < StandardError
end
