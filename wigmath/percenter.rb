#!/usr/bin/env ruby1.9

# == Synopsis 
#   Convert absolute read numbers to read percentages 
#   Divides each value by total
#
# == Examples
#   This command processes seqData.wig:
#      percenter.rb -i seqData.wig -o seqData.percent.wig
#
#   For help use: percenter.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input Wig file to percentize
#	  -t, --total			    Number to divide each value by (optional)	
#   -o, --output        Output percentized Wig file
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
require 'stats'
require 'pickled_optparse'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.wig -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input Wig file" ) { |f| options[:input] = f }
  opts.on( '-t', '--total NUM', "Number to divide each value by (optional)" ) { |n| options[:total] = n.to_f }
  opts.on( '-o', '--output FILE', :required, "Output Wig file (percents)" ) { |f| options[:output] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Initialize Wig file to percentize
wig = WigFile.new(options[:input])

# Number to normalize (divide)
sum = if not options[:total].nil? and options[:total] > 0
	options[:total]
else
	wig.total.to_f
end

# Iterate over the Wig file chromosome-by-chromosome
File.open(options[:output],'w') do |f|
  name = "Sum #{File.basename(options[:output])}"
  desc = "Sum #{File.basename(options[:output])}"
  f.puts Wig.track_header(name, desc) 
  
  wig.each do |chr,values|
    for bp in 0...values.length
      values[bp] /= sum
    end
    
    f.puts Wig.fixed_step(chr, values)
    f.puts values
  end
end