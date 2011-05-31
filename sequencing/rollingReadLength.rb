#!/usr/bin/env ruby1.9

# == Synopsis 
#   Computes the rolling read length across the genome
#
# == Usage 
#   Take BAM reads and computes the avg read length
#   over all bases in the genome
#
#   rollingReadLength.rb -i bowtie.lengths.bed -o readLengths.wig
#
#   For help use: rollingReadLength.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (Bed)
#   -g, --genome        Genome assembly to use
#   -o, --output        Output file with read length averages (Wig)
#
# == Author
#   Timothy Palpant
#
# == Copyright
#   Copyright (c) 2011 Timothy Palpant. Licensed under the GPL v3:
#   http://www.gnu.org/licenses/gpl-3.0.txt

COMMON_DIR = File.expand_path(File.dirname(__FILE__) + '/../../common')
$LOAD_PATH << COMMON_DIR unless $LOAD_PATH.include?(COMMON_DIR)
require 'bundler/setup'
require 'pickled_optparse'
require 'assembly'
require 'wig'
require 'sam'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.bed -o readlengths.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file with reads (Bed)" ) { |f| options[:input] = f }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome ASSEMBLY', "Assembly to use (default: sacCer2)" ) { |g| options[:genome] = g }
  options[:step] = 100_000
  opts.on( '-s', '--step N', "Initial step size to use in base pairs (default: 100,000)" ) { |n| options[:step] = n.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
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


# Initialize a new wig file to hold the output
genome = Assembly.load(options[:genome])
sum = SingleBPData.for_assembly(a)
count = SingleBPData.for_assembly(a)

# Iterate over the reads and map their lengths to the wig file
SAMFile.foreach_read(options[:input]) do |read|
  for bp in read.start-1..read.stop-1
    sum[read.chr][bp] += length
    count[read.chr][bp] += 1
	end
end

# Compute the average length for each base pair
sum.chromosomes.each do |chr|
  for bp in 0...sum[chr].length
    sum[chr][bp] /= count[chr][bp
  end
end

# Set the Wig file metadata
sum.name = "Read Lengths for #{File.basename(options[:input]}"
sum.description = "Read Lengths for #{File.basename(options[:input]}"

# Write to disk
sum.to_wig(options[:output])