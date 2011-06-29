#!/usr/bin/env ruby1.9

# == Synopsis 
#   Compute sequence preferences for paired-end nucleosome reads
#
# == Usage 
#   Take paired-end BAM reads and compute the sequence frequencies up to dinucleotides:
#
#   sequencePreferences.rb -i bowtie.bam -n 2 -o output.txt
#
#   For help use: sequencePreferences.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (BAM)
#   -n, --order         Order (1 = single nucleotide, 2 = dinucleotide, etc)
#   -s, --sequences     Directory with genomic sequences in 2bit format
#   -g, --genome        Genome assembly to use
#   -o, --output        Output file with nucleotide frequencies
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
require 'sam'
require 'genome'
require 'stats'

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
  opts.on( '-s', '--sequences DIR', :required, "Directory with sequences in 2bit format" ) { |d| options[:seqdir] = File.expand_path(d) }
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

# Initialize the Genome
genome_2bit = options[:seqdir] + '/' + options[:genome] + '.2bit'
if not File.exist?(genome_2bit)
  raise "Cannot find directory of FASTA sequences for the genome assembly #{options[:genome]}!"
end 
genome = Genome.new(genome_2bit)

# Initialize the process manager
pm = Parallel::ForkManager.new(options[:threads], {'tempdir' => '/tmp'})

# Callback to get the results from each subprocess
num_bins = 150
counts = Array.new(num_bins, 0)
pm.run_on_finish do |pid, exit_code, ident, exit_signal, core_dump, data|
  # Add the individual chromosome's histogram data to the totals
  for i in 0...num_bins
    counts[i] += data[i]
  end
end


# Process each chromosome in parallel
# Each chromosome in a different parallel process
BAMFile.open(options[:input]) do |bam|
  genome.each do |chr, chr_length|
    # Run in parallel processes managed by ForkManager
    pm.start(chr) and next
    
    puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
    
    chr_hist = Array.new(num_bins, 0)
    
    # Iterate over the read centers on this chromosome, and bin the read length
    bam.each(chr) do |read|
      # Only count paired reads once (use the forward read)
      next if read.paired? and read.crick?
      
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
