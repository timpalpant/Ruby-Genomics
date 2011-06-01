#!/usr/bin/env ruby1.9

# == Synopsis 
#   Find the difference between 2 Wig files
#
# == Usage 
#   Subtract file2.wig from file1.wig:
#
#   difference.rb -m file1.wig -s file2.wig -o output.wig
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
require 'parallelizer'
require 'pickled_optparse'
require 'wig'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -m file1.wig -s file2.wig -p -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-m', '--minuend FILE', :required, "File 1" ) { |f| options[:minuend] = f }
  opts.on( '-s', '--subtrahend FILE', :required, "File 2" ) { |f| options[:subtrahend] = f }
  options[:step] = 200_000
  opts.on( '-c', '--step N', "Chunk size to use in base pairs (default: 200,000)" ) { |n| options[:step] = n.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
      
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end

# Initialize the Wig files to subtract
minuend = WigFile.new(options[:minuend])
subtrahend = WigFile.new(options[:subtrahend])
  
# Validate that both files have the same chromosomes
puts "Validating compatibility"
minuend.chromosomes.each do |chr_id|
  raise "Files have different chromosomes!" unless subtrahend.include?(chr_id)
  raise "Chromosome #{chr_id} has a different length" unless minuend.chr_length(chr_id) == subtrahend.chr_length(chr_id)
end

# Initialize the parallel computation manager
parallelizer = WigComputationParallelizer.new(options[:output], options[:step], options[:threads])

# Run the subtraction on all chromosomes in parallel
parallelizer.run(minuend) do |chr, chunk_start, chunk_stop|
  m_chunk = minuend.query(chr, chunk_start, chunk_stop)
  s_chunk = subtrahend.query(chr, chunk_start, chunk_stop)
  difference = Array.new(m_chunk.length)
  for i in 0...m_chunk.length
    difference[i] = m_chunk[i] - s_chunk[i]
  end
  
  # Return the difference for this chunk
  difference
end
