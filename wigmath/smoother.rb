#!/usr/bin/env ruby1.9

# == Synopsis 
#   Gaussian smooths Wig Files
#
# == Usage 
#   Smooth file1.wig:
#
#   smoother.rb -i file1.wig -o file1.smoothed.wig
#
#   For help use: smoother.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file
#   -s, --stdev         Standard deviation of the Gaussian
#   -w, --window        Window size
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
require 'wig'
require 'stats'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i file1.wig -o file1.smoothed.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file" ) { |f| options[:input] = f }
  options[:sdev] = 10
  opts.on( '-s', '--sdev NUM', "Standard deviation of the Gaussian in base pairs (default 10)" ) { |num| options[:sdev] = num.to_i }
  options[:window_size] = 3
  opts.on( '-w', '--window NUM', "Number of standard deviations +/- to make a window (default 3)" ) { |num| options[:window_size] = num.to_i }
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

  
# Initialize the wig files to smooth
wig = WigFile.new(options[:input])

# Iterate over the Wig file chromosome-by-chromosome
File.open(options[:output],'w') do |f|
  name = "Smoothed #{File.basename(options[:input])}"
  desc = "Smoothed #{File.basename(options[:input])}"
  f.puts Wig.track_header(name,desc)
   
  wig.chromosomes.each do |chr_id|
    chr = wig.chr(chr_id)
    chr.data = chr.gaussian_smooth(options[:sdev],options[:window_size])
    f.puts Wig.fixed_step(chr_id, chr)
    f.puts chr
  end
end