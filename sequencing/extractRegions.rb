#!/usr/bin/env ruby1.9

# == Synopsis 
#   Extract windows from a genome-wide data set and dump to a flat file for 
#   easy import into R or Matlab
#
# == Usage 
#   Dump a list of ORF windows in readcount.marshal to OrfWindows.txt
#
#   extractRegions.rb -i readcount.wig -w windows.bed -o WindowAverage.txt
#
#   For help use: extractRegions.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file to average values (Wig)
#   -w, --windows       List of windows to extract (in Bed format: chrXII  10345  10600)
#   -o, --output        Output file (flat list)
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
require 'bed'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i readcount.wig -w orfs.bed -o ORFWindows.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file to average values (Wig)" ) { |f| options[:input] = f }
  opts.on( '-w', '--window FILE', :required, "List of windows to extract (in Bed format)" ) { |f| options[:windows] = f }
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


puts 'Loading the list of windows'
windows = Bed.load(options[:windows])

# Load the Wig file
wig = Wig.load(options[:input])

puts 'Writing window values to file'
File.open(options[:output],'w') do |f|
  windows.each do |chr,spots|
  	puts "Processing chromosome #{chr}"
  	spots.each do |spot|
  		f.puts wig[chr][spot.start-1..spot.stop-1].to_a.join("\n")
		end
  end
end