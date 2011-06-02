#!/usr/bin/env ruby1.9

# == Synopsis 
#   Log-transform Wig files
#
# == Examples
#   This command processes seqData.wig:
#     logger.rb -i seqData.wig -b 2 -o seqData.log2.wig
#
#   For help use: zscorer.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input Wig file to Z-score
#		-b, --base					Logarithm base (default: 2)
#   -o, --output        Output Z-scored Wig file
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
require 'parallelizer'
require 'wig'
require 'pickled_optparse'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.wig -b 2 -o output.log2.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input Wig file" ) { |f| options[:input] = f }
	options[:base] = 2
	opts.on( '-b', '--base N', "Logarithm base (default: 2)" ) { |n| options[:base] = n.to_i }
  options[:step] = 500_000
  opts.on( '-c', '--step N', "Chunk size to use in base pairs (default: 500,000)" ) { |n| options[:step] = n.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-o', '--output FILE', "Output Wig file (log-transformed)" ) { |f| options[:output] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
  
  # Construct default output filename if not specified
  if options[:output].nil?
    options[:output] = options[:input].split('.')[0..-2].join('.') + ".log#{options[:base]}.wig"
  end
end


# Initialize the Wig file
wig = WigFile.new(options[:input])

# Initialize the parallel computation manager
parallelizer = WigComputationParallelizer.new(options[:output], options[:step], options[:threads])

# Run the subtraction on all chromosomes in parallel
parallelizer.run(wig) do |chr, chunk_start, chunk_stop|
  chunk = wig.query(chr, chunk_start, chunk_stop)
  chunk.map { |value| Math.log(value, options[:base]) }
end