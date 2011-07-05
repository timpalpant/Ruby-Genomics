#!/usr/bin/env ruby1.9

# == Synopsis 
#   Compute sequence preferences for paired-end nucleosome reads
#
# == Usage 
#   Take paired-end input reads and compute the sequence frequencies up to dinucleotides:
#
#   sequencePreferences.rb -i bowtie.input -n 2 -o output.txt
#
#   For help use: sequencePreferences.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (input)
#   -n, --order         Order (1 = single nucleotide, 2 = dinucleotide, etc)
#   -t, --twobit        2bit file with genomic sequences
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
require 'bio-genomic-file'
require 'bio/stats'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.input -o dyads.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input file with reads (input)" ) { |f| options[:input] = f }
  opts.on( '-n', '--order N', :required, "Order of frequencies to compute" ) { |n| options[:order] = n.to_i }
  opts.on( '-t', '--twobit FILE', "Twobit file with genomic reference sequence" ) { |f| options[:twobit] = f }
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

CHUNK_SIZE = 200_000

# Initialize the Genome
genome = Genome.new(options[:twobit])

# What we're searching for
search = ['a', 't', 'c', 'g'].repeated_permutation(options[:order]).map { |p| Bio::Sequence::NA.new(p.join) }

# Initialize the process manager
pm = Parallel::ForkManager.new(options[:threads], {'tempdir' => '/tmp'})

# Callback to get the results from each subprocess
num_bins = 151
half_bins = num_bins / 2
counts = Hash.new
search.each { |n| counts[n] = Array.new(num_bins, 0) }
pm.run_on_finish do |pid, exit_code, ident, exit_signal, core_dump, data|
  # Add the individual chromosome's histogram data to the totals
  data.each do |n, freq|
    for i in 0...num_bins
      counts[n][i] += freq[i]
    end
  end
end


# Process each chromosome in parallel
EntryFile.autodetect(options[:input]) do |input|
  input.chromosomes.each do |chr|
    # Run in parallel processes managed by ForkManager
    pm.start(chr) and next
    
    puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
    
    hist = Hash.new
    search.each { |n| hist[n] = Array.new(num_bins, 0) }
    
    # Search by chunks so that we can cache the sequence
    chunk_start = 1
    while chunk_start <= genome[chr]
      chunk_stop = [chunk_start + CHUNK_SIZE - 1, genome[chr]].min
      puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
      
      # Get the genomic sequence for this chunk
      seq = genome.sequence(chr, chunk_start, chunk_stop)
      
      # Iterate over the reads on this chunk, and tally nucleotide frequencies
      input.each(chr, chunk_start, chunk_stop) do |read|        
        # Get the sequence for this read
        left = read.center - half_bins - chunk_start
        right = read.center + half_bins - chunk_start
        next if left < 0 or right >= CHUNK_SIZE
        read_seq = seq[left..right]
      
        i = options[:order] / 2
        read_seq.window_search(options[:order]) do |subseq|
          hist[subseq][i] += 1
          i += 1
        end
      end
      
      chunk_start = chunk_stop + 1
    end

    # Return the histogram from this subprocess
    pm.finish(0, hist)
  end
end

# Wait for all of the child processes (each chromosome) to complete
pm.wait_all_children

# Write the histogram to the output file
File.open(options[:output], 'w') do |f|
  f.puts "Position\t" + counts.keys.join("\t")
  
  for i in 0...num_bins
    f.puts "#{i-half_bins}\t" + counts.map { |k,v| v[i] }.join("\t")
  end
end
