#!/usr/bin/env ruby1.9

# == Synopsis 
#   Replacing chromosome numbers chr(1-100) with roman numerals
#
# == Usage
#   This command processes genomedata.bed:
#      romanNumeralize.rb -i genomedata.bed -o genomedata.sacCer2.bed
#
#   For help use: romanNumeralize.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file (sacCer1)
#   -o, --output        Output file (sacCer2)
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
require 'roman_numerals'

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
  opts.on( '-i', '--input FILE', :required, "Input file (Arabic integer chromosomes)" ) { |f| options[:input] = f }
  opts.on( '-o', '--output FILE', "Output file (Roman numeral chromosomes)" ) { |f| options[:output] = f }
  
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
  
	# Construct a default output argument if it was not specified
  if options[:output].nil?
  	filename = options[:input].split('.')
  	options[:output] = filename[0..-2].join('.') + '.sacCer2.' + filename[-1] 
  end
end


File.open(options[:output], 'w') do |f|
  File.foreach(options[:input]) do |line|
    # Replace all instances of chr[1-16] with chr[I-XVI]
    converted = line.gsub(/chr[\d]{1,2}/i) do |match|
      # Convert the chromosome number (an integer) to a Roman Numeral
      'chr' + match[3..-1].to_i.to_roman
    end

    f.puts converted
  end
end
