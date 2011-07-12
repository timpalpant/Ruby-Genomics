#!/usr/bin/env ruby1.9

# == Synopsis 
#   Wrapper for abi-dump from sra-toolkit 2.1
#
# == Usage 
#   Dump the reads from the archive test.sra to test.csfasta and test.quals
#
#   abiDumpWrapper.rb -i test.sra -o test.csfasta
#
#   For help use: abiDumpWrapper.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input SRA archive
#   -o, --output        Output CSFASTA/QUALs reads file 1
#   -n, --id            ID for additional output files
#   -d, --directory     Directory for additional output files
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
require 'utils/unix'

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
  opts.on( '-o', '--output FILE', :required, "Output CSFASTA/QUALs" ) { |f| options[:output] = f }
  opts.on( '-n', '--id N', :required, "ID for additional output files" ) { |n| options[:id] = n }
  opts.on( '-d', '--directory DIR', :required, "Directory for additional output files" ) { |d| options[:directory] = d }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Call fastq-dump, output base calls (not colorspace)
output = %x[ abi-dump -F -O #{options[:directory]} #{options[:input]} ]

# Write the output to the Galaxy summary (e.g. "Written 1293 spots")
puts output

# Reorganize the output(s) of abi-dump to fit Galaxy
# See: http://wiki.g2.bx.psu.edu/Admin/Tools/Multiple%20Output%20Files
# abi-dump only lets you specify the folder
output_files = Dir.glob(options[:directory]+'/'+File.basename(options[:input])+'*')
puts "#{output_files.length} output files"
output_files.each_with_index do |f,i|
  puts "Found file #{File.basename(f)}" if ENV['DEBUG']
  extension = File.extname(f)
  if i == 0
    FileUtils.move(f, options[:output])   
  else
    FileUtils.move(f, "#{options[:directory]}/primary_#{options[:id]}_output#{i+1}_visible_#{extension}")
  end
end
