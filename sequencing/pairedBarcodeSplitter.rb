#!/usr/bin/env ruby1.9

# == Synopsis 
#   Split a Fastq file based on paired-end barcodes
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file (mate-paired Fastq)
#   -b, --barcodes      Barcodes file (tab-delimeted)
#   -o, --output        Output file 1
#   -n, --id            Output file 1 ID
#   -d, --directory     Output directory
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
require 'bio'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.fastq -b barcodes.txt -o output1.fastq -i 12323 -d output-dir"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file with mate-paired reads (Fastq)" ) { |f| options[:input] = f }
  opts.on( '-b', '--barcodes FILE', :required, "Barcodes file (tab-delimited)" ) { |f| options[:barcodes] = f }
  opts.on( '-o', '--output FILE', :required, "Output file 1" ) { |f| options[:output] = f }
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


# Load the barcodes
barcodes = Hash.new
File.foreach(options[:barcodes]) do |line|
  next if line.start_with?('#')
  
  entry = line.split("\t")
  raise "Invalid entry in barcodes file: #{line}" if entry.length != 2
  barcodes[entry[0]] = Sequence::NA.new(entry[1])
end

barcode_lengths = barcodes.map { |id,seq| seq.length }.uniq
raise "Barcodes have differening lengths!" if barcode_lengths.length != 1
barcode_length = barcode_lengths.first

# Initialize an output file for each barcode
outputs = Hash.new
i = 1
barcodes.each do |id,seq|
  if i == 1
    outputs[seq] = File.open(options[:output], 'w')
  else
    outputs[seq] = File.open("#{options[:directory]}/primary_#{options[:id]}_output#{i}_visible_fastqsanger", 'w')
  end
  
  i += 1
end

# Iterate over the mate-paired input file and split the reads based on their barcodes
Bio::FlatFile.open(Bio::Fastq, options[:input]) do |f|
  while (read = f.next_entry) and (barcode = f.next_entry)
    str = "@#{read.definition}\n#{read.sequence_string}\n+#{read.definition}\n#{read.quality_string}"
    outputs[barcode.naseq[0, barcode_length]].puts str
  end
end

# Close the output files
outputs.each { |seq,f| f.close }