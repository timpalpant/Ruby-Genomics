#!/usr/bin/env ruby1.9

# == Synopsis 
#   Converts Affymetrix/CisGenome files to BedGraph files
#
# == Usage 
#   affyToBedGraph.rb -i file.txt -o file.bedGraph
#
#   For help use: affyToBedGraph.rb -h
#
# == Options
#   -h, --help          Displays help message
#		-i, --input					Input Affymetrix file
#		-o, --output				Output BedGraph
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
require 'affy'
require 'bedgraph'
require 'pickled_optparse'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
	# Banner at the top of the help screen
	opts.banner = "Usage: #{__FILE__} -i input.txt -o output.bedGraph"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
	opts.on( '-i', '--input FILE', :required, "Input Affymetrix file" ) { |f| options[:input] = f }
	opts.on( '-o', '--output FILE', :required, "Output BedGraph" ) { |f| options[:output] = f }
		
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end

# Initialize the Affymetrix data file
affy = AffyFile.new(options[:input])

# Write to BedGraph format
affy.to_bedgraph(options[:output])