#!/usr/bin/env ruby1.9

# == Synopsis 
#   Wrapper for fasta-get-markov from MEME
#
# == Usage 
#   Generate background Markov model for sequences in seq.fasta:
#
#   memeMarkovWrapper.rb -i seq.fasta -o model.txt
#
#   For help use: memeMarkovWrapper.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input					Input sequences
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

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.fasta -o output.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
	opts.on( '-i', '--input FILE', :required, "Input sequences (FASTA)" ) { |f| options[:input] = f }
	options[:order] = 1
	opts.on( '-m', '--order N', "Markov model order (default = 1)" ) { |n| options[:order] = n.to_i }
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


# Call fasta-get-sequences and swallow the output
# Silence warnings
$VERBOSE = nil
output = %x[ fasta-get-markov -m #{options[:order]} < "\"#{File.expand_path(options[:input])}\"" 2>&1 ].split("\n")

# Hackishly parse the output
markov_line_start = 0
output.each_with_index do |line,line_num|
	# Find the first line that starts with "# order 0"
	# Keep the output past that line
	if line.include?('# order')
		markov_line_start = line_num+1
		break
	end
end

# Write the good output to disk
File.open(options[:output], 'w') do |f|
	f.puts "# order 0"
	f.puts output[markov_line_start..-1].join("\n")
end