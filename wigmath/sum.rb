#!/usr/bin/env ruby1.9

# == Synopsis 
#   Sums Wig Files
#
# == Usage 
#   Sum file1.wig and file2.wig:
#
#   sum.rb file1.wig file2.wig -o wig
#
#   For help use: averager.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -g, --genome        Genome assembly to use (in common/genomes/*)
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
require 'pickled_optparse'
require 'wig'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} file1.wig file2.wig -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
      
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches? or ARGV.length < 2
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end


# Initialize the wig files to add
wigs = ARGV.map { |filename| WigFile.new(filename) }

# Iterate over the Wig files chromosome-by-chromosome
File.open(options[:output],'w') do |f|
  name = "Sum #{File.basename(options[:output])}"
  desc = "Sum #{File.basename(options[:output])}"
  f.puts Wig.track_header(name, desc)  
  
  # TODO: Validate that all the Wigs to average are compatible
  wigs.first.chromosomes.each do |chr|
    sum = wigs.first[chr]
    wigs[1..-1].each do |wig|
      data = wig[chr]
      for bp in 0...sum.length
        sum[bp] += data[bp]
      end
    end
    
    f.puts Wig.fixed_step(chr, sum)
    f.puts sum
  end
end