#!/usr/bin/env ruby1.9

# == Synopsis 
#   Maps the density of read centers
#
# == Usage 
#   Take paired-end Bed reads and count the number of read centers
#   at each loci.
#
#   mapHumanDyads.rb -i bowtie.bed -g hg19 -o dyads.wig
#
#   For help use: mapDyads.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (BAM)
#   -o, --output        Output file with read center density (Wig)
#   -g, --genome        Genome assembly to use (in common/assemblies/*)
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
require 'reference_assembly'
require 'narray'

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
  opts.on( '-g', '--genome NAME', :required, "Genome assembly (default: sacCer2)" ) { |name| options[:genome] = name }
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


# Load the genome assembly
assembly = ReferenceAssembly.load(options[:genome])

File.open(File.expand_path(options[:output]), 'w') do |f|
  name = "Mapped Dyads #{File.basename(options[:input])}"
  f.puts "track type=wiggle_0 name=#{name} description=#{name} visibility=full"

  current_chr = nil
  dyads = nil
  count = 0
  File.foreach(File.expand_path(options[:input])) do |line|
    entry = line.chomp.split("\t")
    chr = entry[0]
    start = entry[1].to_i
    stop = entry[2].to_i
    center = (start + stop) / 2

    if current_chr.nil? or chr != current_chr
      if dyads
        f.puts "fixedStep chrom=#{chr} start=1 step=1 span=1"
        dyads.each { |base| f.puts base }
      end
      current_chr = chr
      puts "Processing chromosome #{chr} (length: #{assembly[chr]})"
      dyads = NArray.sint(assembly[chr])
    end

    dyads[center-1] += 1
    count += 1
    print '.' if count % 1_000_000 == 0
  end
end

# Write the last chromosome
f.puts "fixedStep chrom=#{chr} start=1 step=1 span=1"
dyads.each { |base| f.puts base }
