#!/usr/bin/env ruby1.9

# == Synopsis 
#   Converts BedGraph files to BigWig files
#
# == Usage 
#   bedGraphToWig.rb -i file.bed -o file.bw
#
#   For help use: bedGraphToWig.rb -h
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
require 'bedgraph'
require 'assembly'
require 'pickled_optparse'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
	# Banner at the top of the help screen
	opts.banner = "Usage: #{__FILE__} -i input.bedGraph -g sacCer2 -o output.bw"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
	opts.on( '-i', '--input FILE', :required, "Input BedGraph" ) { |f| options[:input] = f }
	options[:genome] = 'sacCer2'
	opts.on( '-g', '--genome ASSEMBLY', "Genome assembly" ) { |g| options[:genome] = g }
	opts.on( '-o', '--output FILE', :required, "Output BigWig" ) { |f| options[:output] = f }
		
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

# Load the BedGraph data
BedGraphFile.open(options[:input]) do |bedgraph|
  # Write the Wiggle format
  bedgraph.to_bigwig(options[:output], a)
end