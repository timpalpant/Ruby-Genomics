#!/usr/bin/env ruby1.9

# == Synopsis 
#   Compute a DNA structural property on an entire genome and construct a Wig file.
#
# == Examples
#   This command runs Orchid on sacCer2:
#      genomicStructure.rb -g sacCer2 -p orchid -o sacCer2-orchid.wig
#
#   For help use: genomicStructure.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -p, --property      The DNA structural property to compute
#   -g, --genome        Genome assembly
#   -o, --output        Orchid values (Wig)
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
require 'single_bp_data'
require 'structure/dna_property'

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
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome ASSEMBLY', "Genome assembly (default: sacCer2)" ) { |a| options[:genome] = a }
  options[:property] = 'orchid'
  opts.on( '-p', '--property P', "DNA property to compute (default: orchid)" ) { |p| options[:property] = p.downcase }
  opts.on( '-o', '--output FILE', :required, "Output file (Wig)" ) { |f| options[:output] = f }
  
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end

# Load the genome assembly
a = Assembly.load(options[:genome])
seq_file = SingleBPData.for_assembly(a)

# Run structural computation on all chromosomes
puts "Computing #{options[:property]} for all chromosomes" if ENV['DEBUG']
genome.each do |chr_id, bases|  
  puts "Processing chromosome #{chr_id}" if ENV['DEBUG']
  begin
    seq_file[chr_id].data = bases.send(options[:property])
  rescue
    puts "Unknown DNA property: #{options[:property]}" if ENV['DEBUG']
    exit
  end
end

# Write to Wig file
puts "Writing to Wig file" if ENV['DEBUG']
seq_file.to_wig(options[:output])