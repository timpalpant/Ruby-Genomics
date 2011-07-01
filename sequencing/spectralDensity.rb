#!/usr/bin/env ruby1.9

# == Synopsis 
#   Compute the spectral density of a wig dataset in overlapping windows of a given size
#
# == Usage 
#   Take sequencing data in nukes.wig and compute the DFT in 1000bp windows
#   with a step size of 100bp
#
#   spectralDensity.rb -i nukes.wig -w 1000 -s 100 -o output.txt
#
#   For help use: spectralDensity.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with nucleosome sequencing data
#   -w, --window        Window size in base pairs
#   -s, --step          Step size in base pairs
#   -p, --padding       Padding (multiplier) to increase resolution
#   -o, --output        Output file with mean spectral density values
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
require 'fft'
require 'gsl'
include GSL

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.sam -o dyads.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input file (wig)" ) { |f| options[:input] = f }
  opts.on( '-w', '--window N', :required, "Window size in base pairs" ) { |n| options[:window] = n.to_i }
  opts.on( '-s', '--step N', :required, "Step size in base pairs" ) { |n| options[:step] = n.to_i }
  options[:padding] = 0
  opts.on( '-p', '--padding N', "Additional padding multiplier (default: 0)" ) { |n| options[:padding] = n.to_i }
  opts.on( '-o', '--output FILE', :required, "Output file (bedgraph-ish)" ) { |f| options[:output] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end


# Initialize the sequencing data
wig = BigWigFile.new(options[:input])

# Store the mean spectral density
total = Vector[(options[:padding]+1)*options[:window]/2]
count = 0

# Iterate over all the windows and compute the power spectrum
wig.chromosomes.each do |chr|
  puts "Processing chromosome #{chr}" if ENV['DEBUG']
  
  bp = 1
  stop = wig.chr_length(chr)
  while bp + options[:window] < stop
    # Demean
    data = wig.query(chr, bp, bp+options[:window]-1).to_gslv
    data -= data.mean
    
    # Pad with zeros
    padded = Vector[(options[:padding]+1)*data.length]
    padded[0..data.length-1] = data
    
    # Compute the power spectrum and normalize
    p = padded.power_spectrum.to_gslv
    total_power = p.sum
    p /= total_power unless total_power == 0
    
    # Add to the total
    total[0..p.length-1] += p
    count += 1
    
    # Move to the next window
    bp += options[:step]
    # Hack to force GC for large genomes
    GC.start if count % 500 == 0
  end
end

# Compute the average power spectrum (Welch's method for spectral density estimation)
puts "Examined #{count} windows" if ENV['DEBUG']
avg = total / count

# Write to file
File.open(options[:output], 'w') do |f|
  f.puts avg.to_a.join("\n")
end
