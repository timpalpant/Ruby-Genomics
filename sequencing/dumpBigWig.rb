#!/usr/bin/env ruby1.9

# == Synopsis 
#   Dump all data from a BigWig file
#
# == Usage 
#   Dump a list of ORF windows in readcount.bw to OrfWindows.txt
#
#   dumpBigWig.rb -i readcount.bw -o data.txt
#
#   For help use: extractRegions.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file to average values (BigWig)
#   -o, --output        Output file (flat list)
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
require 'forkmanager'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i readcount.wig -w orfs.bed -o ORFWindows.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file (BigWig)" ) { |f| options[:input] = f }
  options[:step] = 500_000
  opts.on( '-s', '--step N', "Step size (default: 500,000bp)" ) { |n| options[:step] = n.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of threads to use (default: 2)" ) { |n| options[:threads] = n.to_i }
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


# Load the Wig file
wig = BigWigFile.new(options[:input])

# Initialize the process manager
pm = Parallel::ForkManager.new(options[:threads])

# Process each chromosome in chunks
# Each chromosome in a different parallel process
wig.chromosomes.each do |chr|
  # Run in parallel processes managed by ForkManager
  pm.start(chr) and next

  puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
  
  File.open(options[:output]+'.'+chr, 'w') do |f|  
    chunk_start = 1
    chr_length = wig.chr_length(chr)
    while chunk_start < chr_length
      # Allocate memory for this chunk
      chunk_stop = [chunk_start+options[:step]-1, chr_length].min
      chunk_size = chunk_stop - chunk_start + 1
      
      # Write this chunk to disk
      f.puts wig.query(chr, chunk_start, chunk_stop).join("\n")
      
      chunk_start = chunk_stop + 1
    end
  end
  
  pm.finish(0)
end


# Wait for all of the child processes (each chromosome) to complete
pm.wait_all_children

# Concatenate all of the individual chromosomes into the output file
tmp_files = wig.chromosomes.map { |chr| options[:output]+'.'+chr }
File.cat(tmp_files, options[:output])

# Delete the individual chromosome files created by each process
tmp_files.each { |filename| File.delete(filename) }
