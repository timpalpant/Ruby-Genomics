#!/usr/bin/env ruby1.9

# == Synopsis 
#   Z-score Wig files
#
# == Examples
#   This command processes seqData.wig:
#     zscorer.rb -i seqData.wig -o seqData.zscored.wig
#
#   For help use: zscorer.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input Wig file to Z-score
#   -o, --output        Output Z-scored Wig file
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
require 'wig'
require 'pickled_optparse'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.wig -o output.zscored.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input Wig file" ) { |f| options[:input] = f }
  opts.on( '-o', '--output FILE', "Output Wig file (Z-scored)" ) { |f| options[:output] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
  
  # Construct default output filename if not specified
  if options[:output].nil?
    options[:output] = options[:input].split('.')[0..-2].join('.') + '.zscored.wig'
  end
end


# Initialize the Wig file
wig = WigFile.new(options[:input])

# Compute mean and standard deviation
mean = wig.mean
puts "Mean: #{mean.to_s(5)}"

stdev = wig.stdev(mean)
puts "StDev: #{stdev.to_s(5)}"
raise "Cannot compute Z-scores for StDev = 0!" if stdev == 0

File.open(options[:input], 'w') do |f|
  name = "Z-Scored #{File.basename(options[:input])}"
  desc = "Z-Scored #{File.basename(options[:input])}"
  f.puts Wig.track_header(name, desc)  
  
  wig.each do |chr,values|
    values.each do |value|
      value = (value-mean) / stdev
    end
    
    f.puts Wig.fixed_step(chr, values)
    f.puts values
  end
end