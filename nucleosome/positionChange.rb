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
#   -w, --wildtype      Wildtype nucleosome calls
#   -m, --mutant        Mutant nucleosome calls
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
require 'stats'
require 'nucleosome'
require 'spot_array'

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
  opts.on( '-w', '--wildtype FILE', :required, "Wildtype nuke calls" ) { |f| options[:wt] = f }
  opts.on( '-m', '--mutant FILE', :required, "Mutant nuke calls" ) { |f| options[:mutant] = f }
  opts.on( '-l', '--loci FILE', :required, "Loci file to use to find +1 nucleosomes (Bed format)" ) { |f| options[:loci] = f }
  options[:num_nukes] = 6
  opts.on( '-n', '--nukes NUM', "Number of nukes to analyze (default: 6)" ) { |n| options[:num_nukes] = n.to_i }
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

if options[:reverse]
  start_nuke = -options[:num_nukes]
  end_nuke = -1
else
  start_nuke = 0
  end_nuke = options[:num_nukes]
end

puts 'Loading Bed list'
loci = Bed.load(options[:loci])
  
puts 'Loading nucleosome calls'
puts '...wildtype'
wildtype = NukeCalls.load(options[:wt])
puts '...mutant'
mutant = NukeCalls.load(options[:mutant])

puts 'Computing change in nucleosome spacing'
position_change = Array.new(options[:num_nukes]) { Array.new }
current_chr, wt_chr_nukes, mutant_chr_nukes = nil, nil, nil
loci.each do |chr,spots|
  puts "Processing chromosome #{chr}"
  
  spots.each do |spot|
    wt_spot_nukes = wildtype[chr].select { |nuke| nuke.start > spot.low and nuke.stop < spot.high }
    mutant_spot_nukes = mutant[chr].select { |nuke| nuke.start > spot.low and nuke.stop < spot.high }
    
    if spot.crick?
      wt_spot_nukes.reverse!
      mutant_spot_nukes.reverse!
    end
    
    for i in start_nuke...end_nuke
      wt_nuke = wt_spot_nukes[i]
      mutant_nuke = mutant_spot_nukes[i]
      next if wt_nuke.nil? or mutant_nuke.nil?
      
      position_change = mutant_nuke.dyad - wt_nuke.dyad
      position_change[i] << position_change unless position_change.abs > 147 # i.e. a nuke call was missed
    end
	end
end

puts 'Writing summary to file'
File.open(options[:output],'w') do |f|
  f.puts "Mean:\t" + position_change.map { |nuke| nuke.mean.to_s(5) }.join("\t")
  f.puts "Median:\t" + position_change.map { |nuke| nuke.median }.join("\t")
  f.puts "LQ:\t" + position_change.map { |nuke| nuke.lower_quartile }.join("\t")
  f.puts "UQ:\t" + position_change.map { |nuke| nuke.upper_quartile }.join("\t")
  f.puts "Min:\t" + position_change.map { |nuke| nuke.min }.join("\t")
  f.puts "Max:\t" + position_change.map { |nuke| nuke.max }.join("\t")
  f.puts "StDev:\t" + position_change.map { |nuke| nuke.stdev.to_s(5) }.join("\t")
  f.puts "Num values:\t" + position_change.map { |nuke| nuke.length }.join("\t")
end

puts 'Writing raw data to file'
File.open(options[:output]+'.raw','w') do |f|
  for i in 0...spots.length
    f.puts position_change.map { |change_arr| change_arr[i] }.map { |nuke| nuke.nil? ? 'NaN' : nuke }.join("\t")
  end
end

global_stop = Time.now
puts "Time elapsed: #{global_stop-global_start}"