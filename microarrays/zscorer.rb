#!/usr/bin/env ruby1.9

# == Synopsis 
#   Z-score GFF files
#
# == Usage 
#   zscorer.rb -i file.gff -o file.zscored.gff
#
#   For help use: zscorer.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input GFF to Z-score
#   -o, --output        Output Z-scored GFF
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
require 'bio-genomic-file'
require 'stats'
require 'pickled_optparse'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  # Banner at the top of the help screen
  opts.banner = "Usage: #{__FILE__} -i input.gff -o output.gff"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input GFF" ) { |f| options[:input] = f }
  opts.on( '-o', '--output FILE', :required, "Output GFF" ) { |f| options[:output] = f }
    
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end


# Load the GFF file
mean, stdev = nil, nil
GFFFile.open(options[:input]) do |gff|
  # Compute mean and standard deviation
  mean = gff.mean
  puts "Mean: #{mean}"

  stdev = gff.stdev
  puts "StDev: #{stdev}"
  raise "Cannot compute Z-scores for StDev = 0!" if stdev == 0
end

# Copy the GFF input file to output, replacing values with Z-scores
print 'Computing Z-scores for each spot...' if ENV['DEBUG']
File.open(options[:output], 'w') do |out|
  File.foreach(options[:input]) do |line|
    # Copy comment lines
    if line.start_with?('#')
      out.write line
    else
      record = line.chomp.split("\t")
      record[5] = (record[5].to_f-avg) / stdev
      out.puts record.join("\t")
    end
  end
end
