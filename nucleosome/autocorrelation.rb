#!/usr/bin/env ruby1.9

# == Synopsis 
#   Compute the autocorrelation of a wig dataset for a list of windows
#
# == Usage 
#   Take sequencing data in nukes.wig and compute the autocorrelation on
#   the list of windows specified in orfs.bed
#
#   autoCorrelation.rb -i nukes.wig -l orfs.bed -o fft.txt
#
#   For help use: autoCorrelation.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with nucleosome sequencing data
#   -l, --loci          Windows to compute the autocorrelation on
#   -o, --output        Output file with autocorrelation values
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
require 'fftw3'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.sam -o dyads.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input file (wig)" ) { |f| options[:input] = f }
  opts.on( '-l', '--loci FILE', :required, "Windows to compute the DFT on (bed)" ) { |name| options[:loci] = name }
  opts.on( '-o', '--output FILE', :required, "Output file (bedgraph-ish)" ) { |f| options[:output] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end


# Initialize the sequencing data
wig = WigFile.autodetect(options[:input])

# Iterate over the list of windows
# and compute the power spectrum for each
File.open(options[:output], 'w') do |f|
  f.puts "ORF\tChromosome\tStart (+1 Nuc)\tStop (3' Nuc)\tL\tAutocorrelation"

  BedFile.foreach(options[:loci]) do |spot|
    next if spot.length <= 1
    
    # Get the data from the Wig
    begin
      data = NArray.to_na(wig.query(spot.chr, spot.low, spot.high).to_a)
    rescue WigError
      next
    end
    
    # Compute the autocorrelation, using the Wiener-Khinchin theorem
    autocorrelation = FFTW3.ifft(FFTW3.fft(data).abs**2).real
    options[:limit] = 200
    #autocorrelation = Array.new(options[:limit])
    #for i in 1..options[:limit]
    #  sum = 0
    #  for j in 0...data.length
    #    sum += data[j] * data[j-i]
    #  end
    #  autocorrelation[i-1] = sum
    #end

    f.puts "#{spot.id}\t#{spot.chr}\t#{spot.start}\t#{spot.stop}\t#{spot.length}\t#{autocorrelation[0...options[:limit]].to_a.join("\t")}"
  end
end
