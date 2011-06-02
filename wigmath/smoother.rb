#!/usr/bin/env ruby1.9

# == Synopsis 
#   Gaussian smooths Wig Files
#
# == Usage 
#   Smooth file1.wig:
#
#   smoother.rb -i file1.wig -o file1.smoothed.wig
#
#   For help use: smoother.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file
#   -s, --stdev         Standard deviation of the Gaussian
#   -w, --window        Window size
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
require 'forkmanager'
require 'pickled_optparse'
require 'fixed_precision'
require 'wig'
require 'stats'
#require 'convolution'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i file1.wig -o file1.smoothed.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file" ) { |f| options[:input] = f }
  options[:sdev] = 20
  opts.on( '-s', '--sdev NUM', "Standard deviation of the Gaussian in base pairs (default 20)" ) { |num| options[:sdev] = num.to_i }
  options[:window_size] = 3
  opts.on( '-w', '--window NUM', "Number of standard deviations +/- to make a window (default 3)" ) { |num| options[:window_size] = num.to_i }
  options[:step] = 200_000
  opts.on( '-c', '--step N', "Chunk size to use in base pairs (default: 200,000)" ) { |n| options[:step] = n.to_i }
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


# Gaussian smoothing requires padding of half_window on either end
padding = options[:sdev] * options[:window_size]

# Set up the Gaussian filter
#g = Filter.gaussian(options[:sdev], options[:window_size])
#padded = Array.new(options[:step], 0)
#for i in 0..g.length
#  padded[i+padded.length/2-g.length/2] = g[i]
#end

# Initialize the wig files to smooth
wig = WigFile.new(options[:input])

# Initialize the process manager
pm = Parallel::ForkManager.new(options[:threads])


# Initialize the parallel computation manager
parallelizer = WigComputationParallelizer.new(options[:output], options[:step], options[:threads])

# Run the subtraction on all chromosomes in parallel
parallelizer.run(wig) do |chr, chunk_start, chunk_stop|
  # Don't pad off the end of the chromosome
  query_start = [1, chunk_start-padding].max
  query_stop = [chunk_stop+padding, wig.chr_length(chr)].min
  
  # Actual padding
  padding_left = chunk_start - query_start
  padding_right = query_stop - chunk_stop
  
  chunk = wig.query(chr, query_start, query_stop)
  smoothed = chunk.gaussian_smooth(options[:sdev], options[:window_size])
  smoothed[padding_left...-padding_right]
end
