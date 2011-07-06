#!/usr/bin/env ruby1.9

# == Synopsis 
#   Convert absolute read numbers to read percentages 
#   Divides each value by total
#
# == Examples
#   This command processes seqData.bw:
#      percenter.rb -i seqData.bw -o seqData.percent.bw
#
#   For help use: percenter.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input BigWig file to percentize
#   -t, --total         Number to divide each value by (optional) 
#   -o, --output        Output percentized BigWig file
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
require 'stats'
require 'pickled_optparse'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.bw -o output.bw"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input BigWig file" ) { |f| options[:input] = f }
  opts.on( '-t', '--total NUM', "Number to divide each value by (optional)" ) { |n| options[:total] = n.to_f }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-g', '--genome ASSEMBLY', :required, "Genome assembly" ) { |g| options[:genome] = g }
  opts.on( '-o', '--output FILE', :required, "Output BigWig file (percent coverages)" ) { |f| options[:output] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Set the number of threads to use
Enumerable.max_threads = options[:threads]

# Initialize Wig file to percentize
wig = WigFile.autodetect(options[:input])

# Number to normalize with (divide by)
sum = if not options[:total].nil? and options[:total] > 0
  options[:total].to_f
else
  wig.total.to_f
end

puts "Normalization factor: #{sum}"

# Initialize the output assembly
assembly = Assembly.load(options[:genome])

# Run the subtraction on all chromosomes in parallel
wig.transform(options[:output], assembly) do |chr, chunk_start, chunk_stop|
  chunk = wig.query(chr, chunk_start, chunk_stop)
  chunk.map { |value| value / sum }
end
