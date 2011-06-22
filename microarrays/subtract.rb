#!/usr/bin/env ruby1.9

# == Synopsis 
#   Subtract one GFF from another and output the result
#
# == Usage 
#   Usage: subtract.rb -x input1.gff -y input2.gff -o difference.gff
#
#   For help use: subtract.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -m, --minuend       Original GFF
#   -s, --subtrahend    Subtract the values of this GFF
#   -o, --output        Output GFF
#
# == Author
#   Timothy Palpant
#
# == Copyright
#   Copyright (c) 2011 Timothy Palpant. Licensed under the GPL v3:
#   http://www.gnu.org/licenses/gpl-3.0.txt

COMMON_DIR = File.expand_path(File.dirname(__FILE__) + '/../../common')
$LOAD_PATH << COMMON_DIR unless $LOAD_PATH.include?(COMMON_DIR)
require 'bundler/setup'
require 'pickled_optparse'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  # Banner at the top of the help screen
  opts.banner = "Usage: subtract.rb -m input.gff -s input2.gff -o difference.gff"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-m', '--minuend FILE', :required, "GFF for original values" ) { |f| options[:m] = f }
  opts.on( '-s', '--subtrahend FILE', :required, "GFF for subtraction values" ) { |f| options[:s] = f }
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


# Read subtraction values in manually so we can hash lookup matches quickly
s = Hash.new
File.foreach(options[:s]) do |line|
  next if line[0] == '#'
  record = line.chomp.split("\t")
  s[record[8]] = record[5].to_f
end

# Copy the GFF input file to output, subtracting values
File.open(options[:output], 'w') do |out|
  File.foreach(options[:m]) do |line|
    # Copy comment lines
    if line[0] == '#'
      out.write line
    else
      record = line.chomp.split("\t")
      next unless s.include?(record[8])
      record[5] = record[5].to_f - s[record[8]]
      out.puts record.join("\t")
    end
  end
end