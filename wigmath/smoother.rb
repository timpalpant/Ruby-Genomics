#!/usr/bin/env ruby1.9

# == Synopsis 
#   Gaussian smooths BigWig Files
#
# == Usage 
#   Smooth file1.bw:
#
#   smoother.rb -i file1.bw -o file1.smoothed.bw
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
require 'pickled_optparse'
require 'bio-genomic-file'
#require 'convolution'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i file1.bw -o file1.smoothed.bw"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file (BigWig)" ) { |f| options[:input] = f }
  options[:sdev] = 20
  opts.on( '-s', '--sdev NUM', "Standard deviation of the Gaussian in base pairs (default 20)" ) { |num| options[:sdev] = num.to_i }
  options[:window_size] = 3
  opts.on( '-w', '--window NUM', "Number of standard deviations +/- to make a window (default 3)" ) { |num| options[:window_size] = num.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-g', '--genome ASSEMBLY', :required, "Genome assembly" ) { |g| options[:genome] = g }
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

# Gaussian smoothing requires padding of half_window on either end
padding = options[:sdev] * options[:window_size]

# Set up the Gaussian filter
#g = Filter.gaussian(options[:sdev], options[:window_size])
#padded = Array.new(options[:step], 0)
#for i in 0..g.length
#  padded[i+padded.length/2-g.length/2] = g[i]
#end

# Initialize the wig files to smooth
wig = WigFile.autodetect(options[:input])

# Initialize the output assembly
assembly = Genomics::Assembly.load(options[:genome])

# Run the subtraction on all chromosomes in parallel
wig.transform(options[:output], assembly, :in_processes => options[:threads]) do |chr, chunk_start, chunk_stop|
  # Don't pad off the end of the chromosome
  query_start = [1, chunk_start-padding].max
  query_stop = [chunk_stop+padding, wig.chr_length(chr)].min
  
  # Actual padding
  padding_left = chunk_start - query_start
  padding_right = query_stop - chunk_stop
  
  chunk = wig.query(chr, query_start, query_stop)
  smoothed = chunk.to_a.gaussian_smooth(options[:sdev], options[:window_size])
  smoothed[padding_left...-padding_right-1]
end
