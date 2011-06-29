#!/usr/bin/env ruby1.9

# == Synopsis 
#   Find the difference between 2 BigWig files
#
# == Usage 
#   Subtract file2.bw from file1.bw:
#
#   difference.rb -m file1.bw -s file2.bw -o output.bw
#
#   For help use: difference.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -m, --minuend       File 1
#   -s, --subtrahend    File 2
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
require 'wig'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -m file1.bw -s file2.bw -p -o output.bw"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-m', '--minuend FILE', :required, "File 1 (BigWig)" ) { |f| options[:minuend] = f }
  opts.on( '-s', '--subtrahend FILE', :required, "File 2 (BigWig)" ) { |f| options[:subtrahend] = f }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-g', '--genome ASSEMBLY', :required, "Genome assembly" ) { |g| options[:genome] = g }
  opts.on( '-o', '--output FILE', :required, "Output file (BigWig)" ) { |f| options[:output] = f }
      
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

# Initialize the BigWig files to subtract
minuend = BigWigFile.new(options[:minuend])
subtrahend = BigWigFile.new(options[:subtrahend])
  
# Validate that both files have the same chromosomes
puts "Validating compatibility"
minuend.chromosomes.each do |chr_id|
  raise "Files have different chromosomes!" unless subtrahend.include?(chr_id)
end

# Initialize the output assembly
assembly = Assembly.load(options[:genome])

# Run the subtraction on all chromosomes in parallel
minuend.transform(options[:output], assembly) do |chr, chunk_start, chunk_stop|
  m_chunk = minuend.query(chr, chunk_start, chunk_stop)
  s_chunk = subtrahend.query(chr, chunk_start, chunk_stop)
  
  # Return the difference for this chunk
  m_chunk - s_chunk
end
