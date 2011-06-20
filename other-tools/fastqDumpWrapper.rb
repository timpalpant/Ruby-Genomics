#!/usr/bin/env ruby1.9

# == Synopsis 
#   Wrapper for fastq-dump from sratoolkit 2.0
#
# == Usage 
#   Dump the reads from the archive test.sra to test.fastq
#
#   fastqDumpWrapper.rb -i test.sra -o test.fastq
#
#   For help use: fastqDumpWrapper.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input					Input SRA archive
#   -o, --output        Output FASTQ reads
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
require 'fileutils'
require 'unix_file_utils'

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
	opts.on( '-i', '--input FILE', :required, "Input SRA archive" ) { |f| options[:input] = f }
  opts.on( '-o', '--output FILE', :required, "Output FASTQ sequences" ) { |f| options[:output] = f }
	
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end

# Get the base directory of the output file
output_dir = File.dirname(File.expand_path(options[:output]))

# Call fastq-dump and swallow the output
# Silence warnings
output = %x[ fastq-dump -SF -O #{output_dir} #{options[:input]} ]

# Write the output (e.g. "Written 1293 spots")
puts output

# Move the output file of fastq-dump to the desired output file name
# since fastq-dump only lets you specify the folder
FileUtils.move(options[:input]+'.fastq', options[:output])
