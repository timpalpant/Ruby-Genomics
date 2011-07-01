#!/usr/bin/env ruby1.9

# == Synopsis 
#   Computes a histogram of read lengths from paired-end sequencing data
#
# == Usage 
#   Take paired-end BAM reads and bin the lengths into bins from 100-200bp
#
#   distogram.rb -i bowtie.bam -r 100:200 -o distogram.txt
#
#   For help use: distogram.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (BAM)
#   -r, --range         Range of bins (bp)
#   -o, --output        Output file with histogram of read lengths
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
require 'forkmanager'
require 'bio-genomic-file'
require 'bio/stats'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.bam -o dyads.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input file with reads (BAM)" ) { |f| options[:input] = f }
  opts.on( '-r', '--range LOW:HIGH', "Range of bins (bp)" ) { |s| options[:range] = s }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome ASSEMBLY', "Assembly to use (default: sacCer2)" ) { |g| options[:genome] = g }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-o', '--output FILE', :required, "Output file (Wig)" ) { |f| options[:output] = f }
  
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

# Load the genome assembly
assembly = Assembly.load(options[:genome])

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
BAMFile.open(options[:input]) do |bam|
  assembly.each do |chr, chr_length|
    # Run in parallel processes managed by ForkManager
    pm.start(chr) and next
    
    puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
    
    chr_hist = Array.new(num_bins, 0)
    
    # Iterate over the read centers on this chromosome, and bin the read length
    bam.each_read(chr) do |read|
      bin = [[read.length, low].max, high].min - low
      chr_hist[bin] += 1
    end

    if ENV['DEBUG']
      puts "Total number of reads on chromosome #{chr}: #{chr_hist.sum}"
      puts "Histogram for chromosome #{chr}:"
      p chr_hist
    end
    
    pm.finish(0, chr_hist)
  end
end

# Wait for all of the child processes (each chromosome) to complete
pm.wait_all_children

# Write the histogram to the output file
File.open(options[:output], 'w') do |f|
  for length in low..high
    f.puts "#{length}\t#{histogram[length-low]}"
  end
end
