#!/usr/bin/env ruby1.9

# == Synopsis 
#   Correlate nucleosome occupancy with occupancy change for overlapping nucleosomes
#
# == Usage 
#   Correlate occupancy vs. change for the nukes in nukes.txt
#
#   occupancyVsChange.rb -n nukes.txt -o nuke-occupancy-change.txt
#
#   For help use: occupancyVsChange.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -n, --nukes			Overlapping nucleosome calls
#   -o, --output        Flat list of occupancy, change in occupancy for all nukes
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
require 'nucleosome'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i nukes.txt -o nuke-distances.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-n', '--nukes FILE', :required, "Overlapping nuke calls" ) { |f| options[:nukes] = f }
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


global_start = Time.now

puts 'Loading overlapping nucleosome calls'
occupancy = Array.new
change = Array.new
File.foreach(options[:nukes]) do |line|
	entry = line.split("\t")
	occ1 = entry[5].to_f
	occ2 = entry[11].to_f
	occ_change = occ2 - occ1
	occupancy << occ1
	change << occ_change
end

puts 'Writing to disk'
File.open(options[:output],'w') do |f|
	for i in 0...occupancy.length
		f.puts "#{occupancy[i]}\t#{change[i]}"
	end
end

global_stop = Time.now
puts "Time elapsed: #{global_stop-global_start}"