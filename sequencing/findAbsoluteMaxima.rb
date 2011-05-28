#!/usr/bin/env ruby1.9

# == Synopsis 
#   Finds the absolute maximum in a list of windows
#
# == Usage 
#  	Finds the peak in the list of ORFS in orfs.bed for the data file data.wig
#
#   findMaxima.rb -i data.wig -l orfs.bed -o orfs-peaks.txt
#
#   For help use: findMaxima.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Data file
#   -l, --loci          Loci to search for peaks in
#   -o, --output        Loci + peak locations
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
require 'stats'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i nukes.txt -o nuke-distances.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Data file" ) { |f| options[:input] = f }
  opts.on( '-l', '--loci FILE', :required, "Loci to search for peaks in (Bed format)" ) { |f| options[:loci] = f }
  opts.on( '-o', '--output FILE', :required, "Output file (Loci + peak locations)" ) { |f| options[:output] = f }
      
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end


# Load the list of windows
loci = Bed.load(options[:loci])
  
# Load the wig data
data = WigFile.new(options[:input])
  
# Find the peak in each locus
loci.each do |chr,spots|
  chr_data = data.chr(chr)

	spots.each do |spot|
    peak_base, peak_value = 0, 0
    values = chr_data.bases(spot.low, spot.high)
		for i in 0...values.length
      if values[i] > peak_value
        peak_base = i
        peak_value = values[i]
      end
    end
    
    spot.value = spot.low + peak_base
	end
end

# Write maxima to file
loci.to_bed(options[:output])