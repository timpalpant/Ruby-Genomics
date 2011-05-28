#!/usr/bin/env ruby1.9

# == Synopsis 
#   Maps the density of read starts
#
# == Usage 
#   Take SAM reads and count the number of read starts
#   at each loci. (Alternative to baseAlignCounts)
#
#   mapStarts.rb -i bowtie.sam -o dyads.wig
#
#   For help use: mapStarts.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (SAM)
#   -o, --output        Output file with read center density (Wig)
#   -g, --genome        Genome assembly to use
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
require 'sam'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.bed -o dyads.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input file with reads (SAM)" ) { |f| options[:input] = f }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome NAME', "Genome assembly (default sacCer2)" ) { |name| options[:genome] = name }
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


a = Assembly.load(options[:genome])
seq_file = SingleBPData.for_assembly(a)

unmapped = 0
SAMFile.foreach_read(options[:input]) do |read|
	begin
		seq_file[read.chr][read.start-1] += 1
	rescue
		unmapped += 1
	end
end

puts "WARN: #{unmapped} unmapped reads" if unmapped > 0

seq_file.to_wig(options[:output])