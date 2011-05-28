#!/usr/bin/env ruby1.9

# == Synopsis 
#   Converts Bed files to Wig files
#
# == Usage 
#   bed2Wig.rb -i file.bed -o file.wig
#
#   For help use: bed2Wig.rb -h
#
# == Options
#   -h, --help          Displays help message
#		-i, --input					Input Bed
#		-g, --genome				Genome assembly
#		-o, --output				Output WIG
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
require 'bed'
require 'assembly'
require 'pickled_optparse'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
	# Banner at the top of the help screen
	opts.banner = "Usage: #{__FILE__} -i input.bed -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
	opts.on( '-i', '--input FILE', :required, "Input Bed" ) { |f| options[:input] = f }
	options[:genome] = 'sacCer2'
	opts.on( '-g', '--genome ASSEMBLY', "Genome assembly" ) { |g| options[:genome] = g }
	opts.on( '-o', '--output FILE', :required, "Output Wig" ) { |f| options[:output] = f }
		
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end

# Load the genome to construct
a = Assembly.load(options[:genome])

# Load the Bed data
bed = Bed.load(options[:input])

# Write the Wiggle format
bed.to_wig(options[:output], a)