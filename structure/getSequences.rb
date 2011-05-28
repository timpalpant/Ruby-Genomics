#!/usr/bin/env ruby1.9

# == Synopsis 
#   Get the sequence underlying a list of loci
#
# == Usage 
#   Get the sequence underlying a list of loci in TSS.txt
#
#   getSequences.rb -l TSS.txt -o orfs-sequences.txt
#
#   For help use: getSequences.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -g, --genome        Genome assembly to use
#   -i, --input         Regions to get the sequences of
#   -o, --output        Sequences (FASTA format)
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

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i ORFS.bed -o sequences.fasta"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome ASSEMBLY', "Assembly to use (default: sacCer2)" ) { |g| options[:genome] = g }
  opts.on( '-i', '--input FILE', :required, "Input file (Bed format)" ) { |f| options[:input] = f }
  opts.on( '-o', '--output FILE', :required, "Output file (FASTA format)" ) { |f| options[:output] = f }
      
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end


global_start = Time.now

puts "Loading genome #{options[:genome]}"
genome = Genome.load(options[:genome])

puts 'Loading list of windows'
loci = Bed.load(options[:input])
  
File.open(options[:output], 'w') do |f|
  loci.each do |chr,spots|
    puts "Processing chromosome #{chr}"
    spots.each do |spot|
      seq = genome[chr].subseq(spot.low, spot.high)
      seq.reverse! if spot.crick?
      f.puts seq.to_fasta
    end
  end
end

global_stop = Time.now
puts "Time elapsed: #{global_stop-global_start}"