#!/usr/bin/env ruby1.9

# == Synopsis 
#   Compute basic statistics on BigWig files
#
# == Examples
#   This command processes seqData.bw:
#     wigstats.rb -i seqData.bw -o seqData.stats.txt
#
#   For help use: wigstats.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input Wig file to Z-score
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
require 'wig'
require 'pickled_optparse'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.bw -o output.stats.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input BigWig file" ) { |f| options[:input] = f }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
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


# Set the number of threads to use
AbstractWigFile.max_threads = options[:threads]

# Initialize the Wig file
wig = BigWigFile.new(options[:input])

# Compute the statistics and write to file
File.open(options[:output], 'w') do |f|
  f.puts wig.to_s
end