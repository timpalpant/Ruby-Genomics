#!/usr/bin/env ruby1.9

# == Synopsis 
#   Finds the first nucleosome from either the 5' or 3' end of a window (such as an ORF)
#
# == Usage 
#  	Finds the first nucleosome for the spots in orfs.bed for the nukes in nukes.txt
#
#   findNuke.rb -i nukes.txt -l orfs.bed -o orfs.nukes.bed
#
#   For help use: findNuke.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Nuke calls
#   -g, --genome        Genome assembly to use
#   -l, --loci          Loci to search for nukes in
#   -o, --output        End-nuke calls
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
require 'bed'
require 'nucleosome'

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
  opts.on( '-i', '--input FILE', :required, "Input file with nuke calls" ) { |f| options[:input] = f }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome ASSEMBLY', "Assembly to use (default: sacCer2)" ) { |g| options[:genome] = g }
  opts.on( '-l', '--loci FILE', :required, "Windows file to find +1 nucleosomes (Bed format)" ) { |f| options[:loci] = f }
  options[:reverse] = false
  opts.on( '-r', '--reverse', "Search from the 3' ends of windows (default: false)" ) { |b| options[:reverse] = b }  
  opts.on( '-o', '--output FILE', :required, "Output file (nuke calls for loci)" ) { |f| options[:output] = f }
	
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end


# Load the windows
loci = Bed.load(options[:loci])

# Mark nucleosomes on the genome
calls = NukeCalls.load(options[:input])

# Pre-sort the nucleosome calls
calls.each { |chr,spots| spots.sort! { |n1,n2| n1.dyad <=> n2.dyad } }

direction = options[:reverse] ? 3 : 5
puts "Finding first nucleosome from #{direction}' end" if ENV['DEBUG']
skipped = 0
invalid_coordinates = 0
loci.each do |chr,spots|
	if not calls.chromosomes.include?(chr)
		skipped += 1
		loci.delete(chr)
		next
	end
	
	spots.each do |spot|
    nukes_in_spot = calls[chr].select { |nuke| nuke.dyad >= spot.low and nuke.dyad <= spot.high }
    if nukes_in_spot.length < 1
			invalid_coordinates += 1
			spots.delete(spot)
			next
		end
		
		# Reverse the nucleosome calls if searching from the 3' end
		nukes_in_spot.reverse! if options[:reverse]
		
    # Take the position of the first nucleosome
    spot.value = nukes_in_spot.first.dyad
	end
end

puts "No nucleosomes for #{skipped} chromosome(s)" if skipped > 0
puts "No nucleosomes for #{invalid_coordinates} coordinate(s)" if invalid_coordinates > 0

# Write to disk
loci.to_bed(options[:output])