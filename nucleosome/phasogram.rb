#!/usr/bin/env ruby1.9

# == Synopsis 
#   Make histogram of distances between dyads
#
# == Usage 
#   Process dyads.bw:
#
#   phasogram.rb -i dyads.bw -r 0:1000 -o dyads.phasogram
#
#   For help use: phasogram.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file (BigWig)
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
require 'forkmanager'
require 'pickled_optparse'
require 'bio-genomic-file'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i dyads.bw -r 0:1000 -o dyads.phasogram"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file (BigWig)" ) { |f| options[:input] = f }
  options[:step] = 500_000
  opts.on( '-c', '--step N', "Chunk size to use in base pairs (default: 500,000)" ) { |n| options[:step] = n.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-r', '--range LOW:HIGH', :required, "Histogram range (bp)" ) { |s| options[:range] = s }
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


# Parse the range
range = options[:range].split(':')
raise "Invalid range given. Range should be of the format LOW:HIGH" if range.length != 2
low = range.first.to_i
high = range.last.to_i
raise "Invalid range given. Range should be of the format LOW:HIGH" if low >= high
num_bins = high - low + 1

# Pad the query ranges so that we get reads with the dyad in the current chunk if it would
# be within the range we're interested in
padding = 2 * high

# Initialize the BigWig dyads file
wig = BigWigFile.new(options[:input])

# Initialize the process manager
pm = Parallel::ForkManager.new(options[:threads], {'tempdir' => '/tmp'})

# Callback to get the results from each subprocess
histogram = Array.new(num_bins, 0)
pm.run_on_finish do |pid, exit_code, ident, exit_signal, core_dump, data|
  # Add the individual chromosome's histogram data to the totals
  for i in 0...num_bins
    histogram[i] += data[i]
  end
end


# Process each chromosome in parallel
# Each chromosome in a different parallel process
wig.chromosomes.each do |chr|
  # Run in parallel processes managed by ForkManager
  pm.start(chr) and next
  
  puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
  
  chr_hist = Array.new(num_bins, 0)
  chunk_start = 1
  chr_length = wig.chr_length(chr)
  while chunk_start < chr_length
    chunk_stop = [chunk_start+options[:step]-1, chr_length].min
    puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']    

    # Don't pad off the end of the chromosome
    query_start = [1, chunk_start-padding].max
    query_stop = [chunk_stop+padding, chr_length].min
    
    # Actual padding
    padding_left = chunk_start - query_start
    padding_right = query_stop - chunk_stop
    
    chunk = wig.query(chr, query_start, query_stop)
    for bp in padding_left...chunk.length-padding_right
      # Get the number of reads a certain distance away
      # Multiply by the number of reads at the current location
      for dist in low..high
        # Left
        chr_hist[dist-low] += chunk[bp] * chunk[bp-dist] if bp-dist >= 0
        # Right
        chr_hist[dist-low] += chunk[bp] * chunk[bp+dist] if bp+dist < chunk.length
      end
    end
    
    chunk_start = chunk_stop + 1
  end
  
  pm.finish(0, chr_hist)
end


# Wait for all of the child processes (each chromosome) to complete
pm.wait_all_children

# Write the histogram to the output file
File.open(options[:output], 'w') do |f|
  for length in low..high
    f.puts "#{length}\t#{histogram[length-low]}"
  end
end
