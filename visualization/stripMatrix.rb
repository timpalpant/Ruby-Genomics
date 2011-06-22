#!/usr/bin/env ruby1.9

# == Synopsis 
#   Strip an aligned matrix for easy import into Matlab
#
# == Usage 
#   Strip the matrix matrix.txt and write to matrix.stripped.txt
#
#   stripMatrix.rb -i matrix.txt -o matrix.stripped.txt
#
#   For help use: stripMatrix.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input matrix
#   -o, --output        Output matrix (stripped)
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
  opts.banner = "Usage: ruby #{__FILE__} -i readcount.wig -l orfs.txt -o orfs-aligned.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input matrix (aligned)" ) { |f| options[:input] = f }
  opts.on( '-o', '--output FILE', :required, "Output matrix (stripped)" ) { |f| options[:output] = f }
	
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end

# Iterate over the lines in the input matrix, stripping column/row headers
line_num = 1
File.open(options[:output],'w') do |out|
	File.foreach(options[:input]) do |line|
		entry = line.chomp.split("\t")
		# Also replace -'s with NaN's for Matlab
		out.puts entry[1..-1].map { |v| (v=='-') ? 'NaN' : v }.join("\t") unless line_num == 1 or entry.first.start_with?('MARKER')
		
		line_num += 1
		puts line_num if line_num % 1_000 == 0 and ENV['DEBUG']
	end
end