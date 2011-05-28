#!/usr/bin/env ruby1.9

# == Synopsis 
#   Mark all relative maxima/minima in a wig file
#		A relative maximum is taken to be a point that is higher than both of its neighbors
#		A relative minimum is taken to be a point that is lower than both of its neighbors
#
# == Usage 
#   Take sequencing data in nukes.wig and make a new wig file with relative minima
#
#   findMinima.rb -i nukes.wig -t minima -o nukes.minima.txt
#
#   For help use: findMinima.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with nucleosome sequencing data
#		-t, --type					minima/maxima
#		-g, --genome				Genome assembly to use (default: sacCer2)
#   -o, --output        Output file with minima marked as 1's and everything else 0
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
require 'assembly'
require 'wig'

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
	options[:type] = 'minima'
	opts.on( '-t', '--type EXTREMA', "Find minima/maxima (default: minima)" ) { |t| options[:type] = t }
	options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome ASSEMBLY', "Assembly to use (default: sacCer2)" ) { |g| options[:genome] = g }
  opts.on( '-o', '--output FILE', :required, "Output file (wig)" ) { |f| options[:output] = f }
	
	# Parse the command-line arguments
	opts.parse!
	options[:type].downcase!
	
	# Validate the required parameters
	if opts.missing_switches? or (options[:type] != 'minima' and options[:type] != 'maxima')
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end

# Boolean flag if we're looking for maxima
maxima = (options[:type] == 'maxima')

# Initialize the sequencing data
wig = WigFile.new(options[:input])

# Initialize a new wig to hold the minima
a = Assembly.load(options[:genome])
extrema = SingleBPData.for_assembly(a)

# Find all extrema
wig.each do |chr,values|
	for bp in 1..values.length-2
		if maxima
			extrema[chr][bp] = 1 if values[bp] > values[bp-1] and values[bp] > values[bp+1]
		else
			extrema[chr][bp] = 1 if values[bp] < values[bp-1] and values[bp] < values[bp+1]
		end
	end
end

# Write to disk
extrema.to_wig(options[:output])