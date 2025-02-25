#!/usr/bin/env ruby1.9

# == Synopsis 
#   Z-score (Big)Wig files
#
# == Examples
#   This command processes seqData.bw:
#     zscorer.rb -i seqData.bw -o seqData.zscored.bw
#
#   For help use: zscorer.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input (Big)Wig file to Z-score
#   -o, --output        Output Z-scored Wig file
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
require 'bio-genomic-file'
require 'pickled_optparse'
require 'reference_assembly'
require 'wig_transform'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.bw -o output.zscored.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input (Big)Wig file" ) { |f| options[:input] = f }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-g', '--genome ASSEMBLY', :required, "Genome assembly" ) { |g| options[:genome] = g }
  opts.on( '-o', '--output FILE', :required, "Output Wig file (Z-scored)" ) { |f| options[:output] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end


# Initialize the Wig file
wig = WigFile.autodetect(options[:input])

# Compute mean and standard deviation
mean = wig.mean
puts "Mean: #{mean.to_s(5)}"

stdev = wig.stdev
puts "StDev: #{stdev.to_s(5)}"
raise "Cannot compute Z-scores for StDev = 0!" if stdev == 0

# Initialize the output assembly
assembly = ReferenceAssembly.load(options[:genome])

# Run the z-scoring on all chromosomes in parallel
wig.transform(options[:output], assembly, :in_processes => options[:threads]) do |chr, chunk_start, chunk_stop|
  chunk = wig.query(chr, chunk_start, chunk_stop)
  chunk.each do |bp,value| 
    chunk.set(bp, (value-mean)/stdev)
  end
  chunk
end
