#!/usr/bin/env ruby1.9

# == Synopsis 
#   Find the average value of windows
#
# == Usage 
#   Find the average value for list of chromosomal windows
#
#   windowStats.rb -i readcount.bw -w windows.bed -o WindowAverage.txt
#
#   For help use: windowStats.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file(s) with genomic data
#   -w, --windows       List of windows to average (in Bed format: chrXII  10345  10600)
#   -s, --statistic     Statistic to compute
#   -o, --output        Output file (flat list)
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
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i readcount.wig -w windows.bed -o WindowAverage.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-w', '--window FILE', :required, "List of windows to average (in Bed format)" ) { |f| options[:windows] = f }
  opts.on( '-s', '--statistic S', "Statistic to compute for each window (default: mean)") { |s| options[:stat] = s }
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
      
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches? or ARGV.length < 1
    puts opts.missing_switches
    puts opts
    exit
  end
  
  # By default, compute the mean
  options[:stat] ||= 'mean'
end

# Load the BigWig files
puts "\nInitializing input data file(s)" if ENV['DEBUG']
inputs = ARGV.map { |input_file| EntryFile.autodetect(input_file) }
puts "Initialized #{inputs.length} files" if ENV['DEBUG']

# Load the windows
puts "Loading window coordinates" if ENV['DEBUG']
loci = Array.new
BedFile.foreach(options[:windows]) { |spot| loci << spot }
puts "Loaded #{loci.length} windows" if ENV['DEBUG']

puts "\nComputing #{options[:stat]} for each window" if ENV['DEBUG']
File.open(options[:output], 'w') do |f|
  basenames = ARGV.map { |input_file| File.basename(input_file) }
  f.puts "#chr\tstart\tstop\tid\t" + basenames.join("\t")

  loci.each do |spot|
    values = inputs.map do |input|
      begin
         value = input.query(spot.chr, spot.low, spot.high).send(options[:stat])
      rescue WigError
         value = 'NaN'
      end

      value
    end

    f.puts "#{spot.to_bed}\t" + values.join("\t")
  end
end

# Close the inputs
inputs.each { |i| i.close }
