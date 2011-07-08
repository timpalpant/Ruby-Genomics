#!/usr/bin/env ruby1.9

# == Synopsis 
#   Maps the density of read centers
#
# == Usage 
#   Take mapDyads.rb Wigs and count the number of read centers
#   at each loci  divided by the read centers in a 147bp window
#   (see: Kaplan et al. 2010)
#
#   conditionalPositioning.rb -i dyads.wig -o conditional.wig
#
#   For help use: conditionalPositioning.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with dyads (Wig)
#   -o, --output        Output file with read center density (Wig)
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
require 'stats'
require 'bio/utils/ucsc'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.bed -o dyads.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file with dyads (Wig)" ) { |f| options[:input] = f }
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


puts "Initializing dyads file" if ENV['DEBUG']
wig = WigFile.autodetect(options[:input])

puts "Computing conditional positioning" if ENV['DEBUG']
File.open(options[:output],'w') do |f|
	name = "#{options[:input]} Conditional Positioning"
  f.write Utils::UCSC::TrackHeader.new(:name => name,
                                       :limits => '0:1e-7')

  wig.each_chr do |chr|
    puts "Processing chromosome #{chr.chr}" if ENV['DEBUG']
    conditional = Genomics::Contig.new(chr.length, chr.chr, 73)
    for bp in 1...chr.length-147
      window = chr.bases(bp, bp+147)
      conditional[bp] = window[73] / window.sum
    end
    
    puts "Writing to disk" if ENV['DEBUG']
    f.puts conditional.to_s
  end
end
