require 'utils/unix'

##
# Tools for calling Perl and Perl scripts from within Ruby code
##
module Perl  
  # Absolute path of the PerlScripts folder
  PERLSCRIPTS = File.expand_path(File.dirname(__FILE__) + '/../resources/perlscripts')
    
  # Run the specified Perl command and return the output
  def self.run(script, args)    
    if not File.exists?(script)
      script = PERLSCRIPTS + '/' + script
      raise "Cannot find Perl script: #{script}" unless File.exists?(script)
    end

    raise "Cannot find Perl interpreter in $PATH" if File.which('perl').nil?
    
    # Execute the Perl script and return the results
    %x[ perl -I#{PERLSCRIPTS} #{script} #{args.join(' ')} ]
  end
end

class PerlError < StandardError
end