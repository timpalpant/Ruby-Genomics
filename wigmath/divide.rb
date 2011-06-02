#!/usr/bin/env ruby1.9

# == Synopsis 
#   Find the division of 2 BigWig files
#
# == Usage 
#   Divide file2.bw into file1.bw:
#
#   divide.rb -m file1.bw -s file2.bw -o output.bw
#
#   For help use: divide.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -1, --dividend      File 1
#   -2, --divisor       File 2
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
  opts.banner = "Usage: ruby #{__FILE__} -1 file1.bw -2 file2.bw -o output.bw"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-1', '--dividend FILE', :required, "File 1" ) { |f| options[:dividend] = f }
  opts.on( '-2', '--divisor FILE', :required, "File 2" ) { |f| options[:divisor] = f }
  options[:step] = 500_000
  opts.on( '-c', '--step N', "Chunk size to use in base pairs (default: 500,000)" ) { |n| options[:step] = n.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-g', '--genome ASSEMBLY', :required, "Genome assembly" ) { |g| options[:genome] = g }
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

# Initialize WigFile files
dividend = BigWigFile.new(options[:dividend])
divisor = BigWigFile.new(options[:divisor])
  
# Validate that both files have the same chromosomes
puts "Validating compatibility"
dividend.chromosomes.each do |chr_id|
  raise "Files have different chromosomes!" unless divisor.include?(chr_id)
  raise "Chromosome #{chr_id} has a different length" unless dividend.chr_length(chr_id) == divisor.chr_length(chr_id)
end

# Initialize the parallel computation manager
parallelizer = BigWigComputationParallelizer.new(options[:output], options[:step], options[:threads])

# Initialize the output assembly
assembly = Assembly.load(options[:genome])

# Run the subtraction on all chromosomes in parallel
parallelizer.run(dividend, assembly) do |chr, chunk_start, chunk_stop|
  dividend_chunk = dividend.query(chr, chunk_start, chunk_stop)
  divisor_chunk = divisor.query(chr, chunk_start, chunk_stop)
  ratio = Array.new(dividend_chunk.length, 0)
  for i in 0...ratio.length
    ratio[i] = dividend_chunk[i] / divisor_chunk[i] unless divisor_chunk[i] == 0
  end
  
  # Return the ratio for this chunk
  ratio
end