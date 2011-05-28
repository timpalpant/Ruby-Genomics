#!/usr/bin/env ruby1.9

# == Synopsis 
#   Compute a DNA structural property on an entire genome and construct a Wig file.
#
# == Examples
#   This command runs Orchid on sacCer2:
#      bindingSiteStructure.rb -g sacCer2 -p orchid -o sacCer2-orchid.wig
#
#   For help use: bindingSiteStructure.rb -h
#
# == Options
#   -h, --help          Displays help message
#		-l, --loci					A list of binding site loci
#   -p, --property      The DNA structural property to compute
#   -g, --genome        Genome assembly
#   -o, --output        Binding site and randomized values
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
require 'genome'
require 'bed'
require 'structure/dna_property'
require 'stats'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -g ASSEMBLY -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
	opts.on( '-l', '--loci FILE', :required, "List of binding sites (Bed)" ) { |f| options[:loci] = f }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome ASSEMBLY', "Genome assembly (default: sacCer2)" ) { |a| options[:genome] = a }
  options[:property] = 'orchid'
  opts.on( '-p', '--property P[P2,P3,...]', Array, "DNA property(s) to compute (default: orchid)" ) { |p| options[:property] = p }
  opts.on( '-o', '--output FILE', :required, "Output file (Wig)" ) { |f| options[:output] = f }
  
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
	
	# Validate the properties
end


puts "Loading genome #{options[:genome]}"
genome = Genome.load(options[:genome])

puts "Loading binding sites"
bs = Bed.load(options[:loci])

# Compute DNA property on all binding sites + random controls
puts "Computing #{options[:property]} for all binding sites (+2 randomizations)"
structure = Array.new
bs.each do |chr,spots|
	print '.'
	spots.each do |spot|
		low = [spot.start, 1].max
		high = [spot.stop, genome[chr].length].min
		
		bs_seq = genome[chr].subseq(low, high)
		random_bs_seq = bs_seq.randomize
		#random_genome_seq = bs_seq.randomize(genome.composition)
		#full_random_seq = bs_seq.randomize({'a'=>1, 't'=>1, 'g'=>1, 'c'=>1})
		#random_locus_seq = 
		
		begin
			bs = bs_seq.send(options[:property]).mean
			random_bs = random_bs_seq.send(options[:property]).mean
			#random_genome = random_genome_seq.send(options[:property]).mean
			#full_random = full_random_seq.send(options[:property]).mean
			#random_locus = 
		rescue NoMethodError
			puts "Unknown DNA property: #{options[:property]}"
			exit
		end
		
		structure << [bs, random_bs]
	end
end
puts 'complete!'

puts "Writing to disk"
File.open(options[:output],'w') do |f|
	f.puts structure.map { |entry| entry.join("\t") }.join("\n")
end