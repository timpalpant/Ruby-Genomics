#!/usr/bin/env ruby1.9

# == Synopsis 
#   Wrapper for WigCorrelate from UCSC-tools
#
# == Usage 
#   Correlate file2.wig with file1.wig:
#
#   wigCorrelate.rb -o output.wig file1.wig file2.wig file3.wig
#
#   For help use: wigCorrelate.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -o, --output        Output file
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
require 'ucsc_tools'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -o output.wig file1.wig file2.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
      
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end


# Call wigCorrelate
output = UCSCTools.wig_correlate(ARGV)

# Parse the output and write to file
File.open(options[:output], 'w') do |f|
  output.split("\n").each do |line|
    entry = line.chomp.split("\t")
    # Clean up the filenames
    file1 = File.basename(entry[0])
    file2 = File.basename(entry[1])
    correlation = entry[2]
    
    # Write to output file
    f.puts "#{file1}\t#{file2}\t#{correlation}"
  end
end