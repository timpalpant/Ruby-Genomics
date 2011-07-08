#!/usr/bin/env ruby1.9

# == Synopsis 
#   Sums BigWig Files
#
# == Usage 
#   Sum file1.bw and file2.bw:
#
#   sum.rb file1.bw file2.bw -o bw
#
#   For help use: averager.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -g, --genome        Genome assembly to use (in common/genomes/*)
#   -o, --output        Output file
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
require 'bio-genomic-file'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} file1.bw file2.bw -o output.bw"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  options[:step] = 500_000
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-g', '--genome ASSEMBLY', :required, "Genome assembly" ) { |g| options[:genome] = g }
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
      
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches? or ARGV.length < 2
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Initialize the wig files to add
wigs = ARGV.map { |filename| WigFile.autodetect(filename) }
# Validate their compatibility
wigs[1..-1].each do |wig|
  wigs.first.chromosomes.each do |chr|
    raise "Wig files do not have the same chromosomes" unless wig.include?(chr)
  end
end

# Initialize the output assembly
assembly = Genomics::Assembly.load(options[:genome])

# Run the subtraction on all chromosomes in parallel
wigs.first.transform(options[:output], assembly, :in_processes => options[:threads]) do |chr, chunk_start, chunk_stop|
  sum = wigs.first.query(chr, chunk_start, chunk_stop)
  wigs[1..-1].each do |wig|
    sum += wig.query(chr, chunk_start, chunk_stop)
  end

  # Return the sum for this chunk
  sum
end
