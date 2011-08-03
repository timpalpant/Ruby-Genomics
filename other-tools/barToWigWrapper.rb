#!/usr/bin/env ruby1.9

# == Synopsis 
#   Wrapper for bar2wig
#
# == Author
#   Timothy Palpant
#
# == Copyright
#   Copyright (c) 2011 Timothy Palpant. Licensed under the GPL v3:
#   http://www.gnu.org/licenses/gpl-3.0.txt

COMMON_DIR = File.expand_path(File.dirname(__FILE__) + '/../common')
$LOAD_PATH << COMMON_DIR unless $LOAD_PATH.include?(COMMON_DIR)
require 'bundler/setup'
require 'pickled_optparse'
require 'fileutils'
require 'python'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.bar -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input BAR file" ) { |f| options[:input] = f }
  opts.on( '-o', '--output FILE', :required, "Output Wig file" ) { |f| options[:output] = f }
  opts.on( '-t', '--threshold F', "Minimum threshold (default: 0.0)" ) { |i| options[:threshold] = i.to_f }
  opts.on( '-s', '--step N', "Step size (default = 1)" ) { |n| options[:step] = n.to_i }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Call bar2wig
bar2wig_script = File.expand_path(File.dirname(__FILE__) + '/bar2wig.py')
args = [options[:input]]
args << '-t' << options[:threshold] if options.include?(:threshold)
args << '-s' << options[:step] if options.include?(:step)
Python.run(bar2wig_script, args)

# Move the output file
default_output = options[:input] + '.wig'
FileUtils.move(default_output, options[:output])
