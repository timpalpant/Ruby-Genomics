#!/usr/bin/env ruby1.9

# == Synopsis 
#   Compute the power spectrum of a wig dataset for a list of windows
#
# == Usage 
#   Take sequencing data in nukes.wig and compute the DFT on
#		the list of windows specified in orfs.bed
#
#   powerSpectrum.rb -i nukes.wig -l orfs.bed -o fft.txt
#
#   For help use: powerSpectrum.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with nucleosome sequencing data
#   -l, --loci					Windows to compute the DFT on
#   -o, --output        Output file with DFT values
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
require 'fft'

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
  opts.on( '-l', '--loci FILE', :required, "Windows to compute the DFT on (bed)" ) { |name| options[:loci] = name }
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


# Load the list of windows
loci = Bed.load(options[:loci])

# Initialize the sequencing data
wig = WigFile.new(options[:input])

# Iterate over the list of windows
# and compute the power spectrum for each
crystal, bistable, other = 0, 0, 0
loci.each do |chr,spots|
	puts "chromosome #{chr}" if ENV['DEBUG']
	# Load the data for the current chromosome from disk
	begin
		chr_data = wig.chr(chr)
	rescue GenomicIndexError
		next
	end
	
	# Process the current chromosome
	spots.each do |spot|
		next if spot.length <= 1
		
		# Compute the power spectrum
		begin
			p = chr_data.bases(spot.start, spot.stop).power_spectrum
		rescue GenomicIndexError
			next
		end
		
		# Normalize
		p /= p.sum
		
		# Decide whether the window is crystal, bistable, or other
		# based on definitions in Vaillant et al. 2010
		sorted_indices = p.sort_index.reverse
		sorted = p[sorted_indices]
		
		num_nukes = sorted_indices.to_v + 1
		period = spot.length * num_nukes**-1
		ratio = sorted[0] / sorted[1]
		mean_period = (period[0]+period[1]) / 2
		last_freq = [100, p.length-1].min
		if sorted[0] > 0.5 and ratio > 3 and period[0] > 125 and period[0] < 210
			crystal += 1
			spot.value = "Crystal"
		elsif ratio > 0.5 and ratio < 2 and mean_period > 147 and mean_period < 180
			bistable += 1
			spot.value = "Bistable"
		else
			other += 1
			spot.value = "Other"
		end
		
		spot.value << "\t#{sorted[0]}\t#{num_nukes[0]}\t#{period[0]}\t#{sorted[1]}\t#{num_nukes[1]}\t#{period[1]}\t\t#{p[0..last_freq].to_a.join("\t")}" 
	end
end

puts "Crystal: #{crystal}"
puts "Bistable: #{bistable}"
puts "Other: #{other}"

# Write the power spectra to disk as a comma-separated list in the 5th column
File.open(options[:output], 'w') do |f|
	f.puts "ORF\tChromosome\tStart (+1 Nuc)\tStop (3' Nuc)\tL\tType\tPS1\tNumNuc1\tPeriod1\tPS2\tNumNuc2\tPeriod2\t\tPS values"
	loci.each do |chr,spots|
		spots.each do |spot|
			f.puts "#{spot.id}\t{chr}\t#{spot.start}\t#{spot.stop}\t#{(spot.start-spot.stop).abs}\t#{spot.value}"
		end
	end
end