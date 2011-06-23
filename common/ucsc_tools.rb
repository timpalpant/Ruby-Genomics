require 'unix_file_utils'
require 'stringio'

##
# Wrap UCSC tools programs for using them in Ruby scripts
##
module UCSCTools
  # Run the specified executable and return the output
  def self.run(command)
    command_line = command.split(' ')
    program = command_line.first
    args = command_line[1..-1].join(' ')

    raise "Cannot find executable: #{program} in $PATH" if File.which(program).nil?
    
    # Execute the program and return the results
    %x[ #{program} #{args} ]
  end
	
  def self.wig_correlate(files)
    run("wigCorrelate #{files.join(' ')}")
  end
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
    
    track = self.new
    pos = 0
    while (equals_pos = line.index('=', pos))
      begin
        # Look back from the equals position until there is a space to get the token key
        cursor = equals_pos - 1
        cursor -= 1 while line[cursor] != ' '
        key = line[cursor+1, equals_pos-1]
        
        # Look forward from the equals position until there is a space to get the token value
        cursor = equals_pos + 1
        cursor += 1 while line[cursor] != ' '
        value = line[equals_pos+1, cursor-1]
        
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
