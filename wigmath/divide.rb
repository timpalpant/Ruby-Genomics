#!/usr/bin/env ruby1.9

# == Synopsis 
#   Find the division of two (Big)Wig files
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
require 'pickled_optparse'
require 'bio-genomic-file'
require 'reference_assembly'
require 'wig_transform'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -1 file1.bw -2 file2.bw -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-1', '--dividend FILE', :required, "File 1" ) { |f| options[:dividend] = f }
  opts.on( '-2', '--divisor FILE', :required, "File 2" ) { |f| options[:divisor] = f }
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
dividend = WigFile.autodetect(options[:dividend])
divisor = WigFile.autodetect(options[:divisor])
  
# Validate that both files have the same chromosomes
puts "Validating compatibility" if ENV['DEBUG']
dividend.chromosomes.each do |chr_id|
  raise "Files have different chromosomes!" unless divisor.include?(chr_id)
end

# Initialize the output assembly
assembly = ReferenceAssembly.load(options[:genome])

# Run the subtraction on all chromosomes in parallel
dividend.transform(options[:output], assembly, :in_processes => options[:threads]) do |chr, chunk_start, chunk_stop|
  dividend_chunk = dividend.query(chr, chunk_start, chunk_stop)
  divisor_chunk = divisor.query(chr, chunk_start, chunk_stop)
  
  # Return the ratio for this chunk
  dividend_chunk / divisor_chunk
end
