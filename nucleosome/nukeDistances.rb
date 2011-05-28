#!/usr/bin/env ruby1.9

# == Synopsis 
#   Calculates internucleosome distances from a nuke calls file
#
# == Usage 
#   Calculate internucleosome distances for the nukes in nukes.txt
#
#   nukeDistances.rb -n nukes.txt -o nuke-distances.txt
#
#   For help use: nukeDistances.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Nuke calls
#   -g, --genome        Genome assembly to use
#   -o, --output        Flat list of nuke distances
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
require 'assembly'
require 'single_bp_data'
require 'spot'
require 'stats'
require 'roman_numerals'

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
  opts.on( '-i', '--input FILE', :required, "Input file with nuke calls" ) { |f| options[:input] = f }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome ASSEMBLY', "Assembly to use (default: sacCer2)" ) { |g| options[:genome] = g }
  opts.on( '-l', '--loci FILE', :required, "Loci file to use to find +1 nucleosomes (Bed format)" ) { |f| options[:loci] = f }
  options[:cumulative] = false
  opts.on( '-c', '--cumulative', "Should distances computed by cumulative (default: false)" ) { |b| options[:cumulative] = b }
  options[:reverse] = false
  opts.on( '-r', '--reverse', "Search from the ends of loci (default: false)" ) { |b| options[:reverse] = b }  
  opts.on( '-o', '--output FILE', :required, "Output file (internucleosome distances)" ) { |f| options[:output] = f }
      
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

puts 'Loading loci list'
loci = BedGraph.load(options[:loci])

puts 'Initializing genome file'
a = Assembly.load(options[:genome])
wig = SingleBPData.for_assembly(a)
  
puts 'Marking all nucleosomes on the genome'
line_num = 0
File.foreach(options[:input]) do |line|
  # Skip header line
  line_num += 1
  next if line_num == 1
  
  entry = line.split("\t")
  chr = entry[0]
  center = entry[5].to_i
  
  # Mark the nuke position
  wig[chr][center] += 1
end

puts 'Computing internucleosomal distances'
distances = Array.new
loci.each do |spot|
  low = [spot.start, spot.stop].min
  high = [spot.start, spot.stop].max
  values = wig[spot.chr-1][low-1..high-1]
  if options[:reverse]
    values.reverse! if spot.stop > spot.start
  else
    values.reverse! if spot.start > spot.stop
  end

  distance = 0
  nuke = 0
  values.each do |value|
    if value > 0
      distances[nuke] ||= Array.new
      distances[nuke] << distance
      nuke += 1
      distance = 0 if nuke == 1 or not options[:cumulative]
    end
    
    distance += 1
  end
end

puts 'Writing summary to file'
File.open(options[:output],'w') do |f|
  f.puts "Mean:\t" + distances.map { |nuke| nuke.mean.to_s(5) }.join("\t")
  f.puts "Median:\t" + distances.map { |nuke| nuke.median }.join("\t")
  f.puts "LQ:\t" + distances.map { |nuke| nuke.lower_quartile }.join("\t")
  f.puts "UQ:\t" + distances.map { |nuke| nuke.upper_quartile }.join("\t")
  f.puts "Min:\t" + distances.map { |nuke| nuke.min }.join("\t")
  f.puts "Max:\t" + distances.map { |nuke| nuke.max }.join("\t")
  f.puts "StDev:\t" + distances.map { |nuke| nuke.stdev.to_s(5) }.join("\t")
  f.puts "Num values:\t" + distances.map { |nuke| nuke.length }.join("\t")
end

global_stop = Time.now
puts "Time elapsed: #{global_stop-global_start}"