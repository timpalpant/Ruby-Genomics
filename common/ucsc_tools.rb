require 'unix_file_utils'
require 'bio'
require 'stringio'

##
# Wrap UCSC tools programs for using them in Ruby scripts
##
module UCSCTools
  def self.wig_correlate(files)
    %x[ wigCorrelate #{files.map { |f| File.expand_path(f) }.join(' ')} ]
  end
  
  def self.bedgraph_to_bigwig(input_file, assembly_file, output_file)
    %x[ bedGraphToBigWig #{File.expand_path(input_sorted)} #{File.expand_path(assembly_file)} #{File.expand_path(output_file)} ]
  end
  
  def self.bed_to_bigbed
    raise "Not yet implemented"
  end
  
  def self.bigbed_info
    raise "Not yet implemented"
  end
  
  def self.bigbed_summary
    raise "Not yet implemented"
  end
  
  def self.bigbed_to_bed
    raise "Not yet implemented"
  end
  
  def self.bigwig_info(f)
    %x[ bigWigInfo -chroms #{File.expand_path(f)} ].split("\n")
  end
  
  def self.bigwig_summary(f, chr, start, stop, num_values, type = 'mean')
    # Data is 0-indexed and half-open
    output = %x[ bigWigSummary -type=#{type} #{File.expand_path(f)} #{chr} #{start-1} #{stop} #{num_values} 2>&1 ]
    raise UCSCToolsError, "BigWig does not contain data for the interval #{chr}:#{start}-#{stop}" if output.start_with?('no data in region')

    values = output.split(' ')
    raise UCSCToolsError, "bigWigSummary did not return the expected number of values!" if values.length != num_values

    return values.map { |v| v.to_f unless v == 'n/a' or v == 'nan' }
  end
  
  def self.bigwig_to_bedgraph(input_file, output_file)
    %x[ bigWigToBedGraph #{File.expand_path(input_file)} #{File.expand_path(output_file)} ]
  end
  
  def self.bigwig_to_wig(input_file, output_file)
    %x[ bigWigToWig #{File.expand_path(input_file)} #{File.expand_path(output_file)} ]
  end
  
  def self.wig_to_bigwig(input_file, assembly_file, output_file)
    %x[ wigToBigWig -clip #{File.expand_path(input_file)} #{File.expand_path(assembly_file)} #{File.expand_path(output_file)} ]
  end
  
  def self.twobit_info(twobit_file)
    result = Hash.new
    
    %x[ twoBitInfo #{File.expand_path(twobit_file)} stdout ].split("\n").each do |line| 
      entry = line.split("\t")
      result[entry.first] = entry.last.to_i
    end
    
    return result
  end
  
  def self.twobit_to_fa(twobit_file, chr, start = nil, stop = nil)
    # TwoBit is 0-indexed and half-open
    query_string = chr
    query_string += ":#{start-1}-#{stop}" if start and stop
    
    s = StringIO.new
    IO.popen("twoBitToFa #{File.expand_path(twobit_file)}:#{query_string} stdout") do |output|
      output.each do |line|
        # Skip the header line
        next if line.start_with?('>')
        s << line.chomp
      end
    end

    seq = Bio::Sequence::NA.new(s.string)
    if (start and stop) and seq.length != (stop-start+1)
      raise UCSCToolsError, "twoBitToFa did not return the expected sequence length!" 
    end

    return seq
  end
end

class UCSCToolsError < StandardError
end

##
# A track header line for the UCSC Genome Browser
# For the spec, see: http://genome.ucsc.edu/goldenPath/help/wiggle.html
##
class UCSCTrackHeader
  attr_accessor :type, :name, :description, :visibility, :color, :alt_color, :priority, :auto_scale, :always_zero, :grid_default, :max_height_pixels, :graph_type, :view_limits, :y_line_mark, :y_line_on_off, :windowing_function, :smoothing_window

  # TODO: Validate the parameters
  def initialize(opts = {})
    @type = opts[:type]
    @name = opts[:name]
    @description = opts[:description]
    @visibility = opts[:visibility]
    @color = opts[:color]
    @alt_color = opts[:alt_color]
    @priority = opts[:priority]
    @auto_scale = opts[:auto_scale]
    @always_zero = opts[:always_zero]
    @grid_default = opts[:grid_default]
    @max_height_pixels = opts[:max_height_pixels]
    @graph_type = opts[:graph_type]
    @view_limits = opts[:view_limits]
    @y_line_mark = opts[:y_line_mark]
    @y_line_on_off = opts[:y_line_on_off]
    @windowing_function = opts[:windowing_function]
    @smoothing_window = opts[:smoothing_window]
  end
  
  def set(key, value)
    puts "Setting UCSC track #{key}:#{value}" if ENV['DEBUG']

    case key
      when 'type'
        @type = value
      when 'name'
        @name = value
      when 'description'
        @description = value
      when 'visibility'
        @visibility = value
      when 'color'
        @color = value
      when 'altColor'
        @alt_color = value
      when 'priority'
        @priority = value
      when 'autoScale'
        @auto_scale = value
      when 'alwaysZero'
        @always_zero = value
      when 'gridDefault'
        @grid_default = value
      when 'maxHeightPixels'
        @max_height_pixels = value
      when 'graphType'
        @graph_type = value
      when 'viewLimits'
        @view_limits = value
      when 'yLineMark'
        @y_line_mark = value
      when 'yLineOnOff'
        @y_line_on_off = value
      when 'windowingFunction'
        @windowing_function = value
      when 'smoothingWindow'
        @smoothing_window = value
      else
        raise UCSCTrackHeaderError, "Unknown UCSC track header key #{key}"
    end
  end
  
  # Parse the tokens in a track line into a UCSCTrackHeader object
  def self.parse(line)
    unless line.chomp.start_with?('track')
      raise UCSCTrackHeaderError, "Not a valid UCSC Genome Browser track line"
    end

    puts "Parsing track line: #{line.chomp!}" if ENV['DEBUG']
    
    track = self.new
    pos = 0
    while (equals_pos = line.index('=', pos))
      begin
        # Look back from the equals position until there is a space to get the token key
        cursor = equals_pos - 1
        stop = cursor
        cursor -= 1 until cursor == 0 or line[cursor] == ' '
        start = cursor + 1
        key = line[start..stop]
        
        # Look forward from the equals position until there is a space to get the token value
        cursor = equals_pos + 1
        quoted = (line[cursor] == '"' or line[cursor] == "'")
        cursor += 1 if quoted
        start = cursor
        cursor += 1 until cursor == line.length or (not quoted and line[cursor] == ' ') or (quoted and (line[cursor] == '"' or line[cursor] == "'"))
        stop = cursor - 1
        value = line[start..stop]        

        # Store the token key-value in the UCSCTrackHeader object
        begin
          track.set(key, value)
        rescue UCSCTrackHeaderError
          puts "Unknown UCSC track header key: #{key}, ignoring" if ENV['DEBUG']
        end
      rescue
        puts "Malformed UCSC track header line" if ENV['DEBUG']
      ensure
        # Move to the next token
        pos = equals_pos + 1
      end
    end
    
    return track
  end

  def to_s  
  	str = StringIO.new
    str << "track"
    
    str << " type=#{@type}" if @type
    str << " name='#{@name}'" if @name
    str << " description='#{@description}'" if @description
    str << " autoScale=#{@auto_scale}" if @auto_scale
    str << " visibility=#{@visibility}" if @visibility
    str << " viewLimits=#{@view_limits}" if @view_limits
    str << " color=#{@color}" if @color
    str << " altColor=#{@alt_color}" if @alt_color
    str << " priority=#{@priority}" if @priority
    str << " alwaysZero=#{@always_zero}" if @always_zero
    str << " gridDefault=#{@grid_default}" if @grid_default
    str << " maxHeightPixels=#{@max_height_pixels}" if @max_height_pixels
    str << " graphType=#{@graph_type}" if @graph_type
    str << " yLineMark=#{@y_line_mark}" if @y_line_mark
    str << " yLineOnOff=#{@y_line_on_off}" if @y_line_on_off
    str << " windowingFunction=#{@windowing_function}" if @windowing_function
    str << " smoothingWindow=#{@smoothing_window}" if @smoothing_window
    
    return str.string
  end
end

class UCSCTrackHeaderError < StandardError
end
