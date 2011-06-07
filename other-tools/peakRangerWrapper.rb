#!/usr/bin/env ruby1.9

# == Synopsis 
#   Wrapper for Peak Ranger
#
# == Usage 
#   Call peaks on ChIP data vs input control (required):
#
#   peakRangerWrapper.rb -s sample.bam -c control.bam -o peaks.txt
#
#   For help use: peakRangerWrapper.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -s, --sample        ChIP (BAM)
#   -c, --control       Control (BAM)
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
require 'fileutils'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i input.fasta -o output.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
	opts.on( '-s', '--sample FILE', :required, "ChIP (BAM)" ) { |f| options[:sample] = f }
	opts.on( '-c', '--control FILE', :required, "Control (BAM)" ) { |f| options[:control] = f }
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

# Call ranger and swallow the output
output = %x[ ranger --format=bam -d #{options[:sample]} -c #{options[:control]} -o #{options[:output]} #{ARGV.join(' ')} ]

# Delete unwanted output files
peaks_bed = options[:output] + '_peaks.bed'
peaks_with_region = options[:output] + '_peaks_with_region.bed'
raw = options[:output] + '_raw'
regions = options[:output] + '_regions'
valley_bed = options[:output] + '_valley.bed'
valleys_with_region = options[:output] + '_valleys_with_region.bed'
[peaks_bed, raw, regions, valley_bed, valleys_with_region].each { |f| File.delete(f) }

# Move the output to keep into its expected location
FileUtils.move(peaks_with_region, options[:output])