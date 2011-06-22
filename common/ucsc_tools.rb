require 'unix_file_utils'

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
	
  # Correlate multiple Wiggle files
	def self.wig_correlate(files)
		run("wigCorrelate #{files.join(' ')}")
	end
  
  # Construct a track line for the UCSC Genome Browser
  # TODO: Validate that the arguments are valid
  def self.track_header(opts = {})
    # By default, if no description is specified and a name is, use the name as the description
    opts[:description] ||= opts[:name] if opts[:name]
  
  	str = StringIO.new("track")
    
    str << " type=#{opts[:type]}" if opts[:type]
    str << " name='#{opts[:name]}'" if opts[:name]
    str << " description='#{opts[:description]}'" if opts[:description]
    str << " autoScale=#{opts[:autoscale]}" if opts[:autoscale]
    str << " visibility=#{visibility}" if opts[:visibility]
    str << " viewLimits=#{opts[:limits]}" if opts[:limits]
    
    return str.string
  end
end