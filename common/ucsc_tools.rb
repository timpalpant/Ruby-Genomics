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
  
  def self.wig_correlate(files)
    run("wigCorrelate #{files.join(' ')}")
  end
end