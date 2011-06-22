#!/usr/bin/env ruby1.9

# == Synopsis 
#   Get the base pair distribution of a set of equal-length sequences
#
# == Usage 
#   Get the base pair distribution of sequences in seq.fasta
#
#   sequenceDistribution.rb -i seq.fasta -o seq-dist.txt
#
#   For help use: sequenceDistribution.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         FASTA file with equal-length sequences
#   -o, --output        Nuke comparison results
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
require 'fixed_precision'
require 'genome'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i seq.fasta -o seq-distribution.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file (Fasta format)" ) { |f| options[:input] = f }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome NAME', "Genome assembly (default sacCer2)" ) { |name| options[:genome] = name }
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


global_start = Time.now

# Load the genome
genome = Genome.load(options[:genome])

# Load the sequences as a (# seq x 4) matrix with rows = base pair
# and columns (Hash) being the base pair count
distribution = Array.new
num_sequences = 0.0
Bio::FlatFile.foreach(Bio::FastaFormat, options[:input]) do |entry|
  num_sequences += 1
  entry.naseq.split(//).each_with_index do |base,i|
    distribution[i] ||= Hash.new(0)
    distribution[i][base] += 1
  end
end

# Compute the frequency of each base pair across the sequence
a = distribution.map { |base| base['a'] / num_sequences }
t = distribution.map { |base| base['t'] / num_sequences }
g = distribution.map { |base| base['g'] / num_sequences }
c = distribution.map { |base| base['c'] / num_sequences }

File.open(options[:output], 'w') do |f|
  f.puts "Genome-wide:\tA\tT\tG\tC"
  f.puts "\t#{genome.a_content.to_f.to_s(5)}\t#{genome.t_content.to_f.to_s(5)}\t#{genome.g_content.to_f.to_s(5)}\t#{genome.c_content.to_f.to_s(5)}"
  f.print "\n"
  
  f.puts "A:\t" + a.map { |elem| elem.to_s(5) }.join("\t")
  f.puts "T:\t" + t.map { |elem| elem.to_s(5) }.join("\t")
  f.puts "G:\t" + g.map { |elem| elem.to_s(5) }.join("\t")
  f.puts "C:\t" + c.map { |elem| elem.to_s(5) }.join("\t")
end

global_stop = Time.now
puts "Time elapsed: #{global_stop-global_start}"