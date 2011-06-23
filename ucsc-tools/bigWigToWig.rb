#!/usr/bin/env ruby1.9

# == Synopsis 
#   Wrapper for bigWigToWig from UCSC-tools
#
# == Usage 
#   Correlate file2.wig with file1.wig:
#
#   bigWigToWig.rb -i input.wig -o output.bw
#
#   For help use: bigWigToWig.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file
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
require 'wig'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.wig -o output.bw"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file (Wig)" ) { |f| options[:input] = f }
  opts.on( '-o', '--output FILE', :required, "Output file (BigWig)" ) { |f| options[:output] = f }
      
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end

BigWigFile.to_wig(options[:input], options[:output])
