require 'utils/unix'

##
# Tools for calling Python and Python scripts from within Ruby code
##
module Python  
  # Absolute path of the PerlScripts folder
  PYTHONSCRIPTS = File.expand_path(File.dirname(__FILE__) + '/../resources/pythonscripts')
    
  # Run the specified Perl command and return the output
  def self.run(script, args)    
    if not File.exists?(script)
      script = PYTHONSCRIPTS + '/' + script
      raise "Cannot find Python script: #{script}" unless File.exists?(script)
    end

    raise "Cannot find Python interpreter in $PATH" if File.which('python').nil?
    
    # Execute the Python script and return the results
    %x[ python -I#{PYTHONSCRIPTS} #{script} #{args.join(' ')} ]
  end
end
