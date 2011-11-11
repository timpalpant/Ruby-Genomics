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
require 'bio-genomic-file'
include Bio

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

CHUNK_SIZE = 1_000_000

# Load the Wig file
wig = WigFile.autodetect(options[:input])

# Process each chromosome in chunks
File.open(options[:output], 'w') do |f|
  wig.chromosomes.each do |chr|
    start = wig.chr_start(chr)
    stop = wig.chr_stop(chr)

    bp = start
    while bp < stop
      query_start = bp
      query_stop = [bp+CHUNK_SIZE-1, stop].min
      data = wig.query(chr, query_start, query_stop)
      data.each do |bp, value|
        f.puts value
      end
    end
  end
end
