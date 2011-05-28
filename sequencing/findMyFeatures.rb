#!/usr/bin/env ruby1.9

# == Synopsis 
#   Take peak calls and determine which feature each peak is in.
#
# == Usage 
#   Find the features associated with Rap1 peaks:
#
#   findMyFeatures.rb -i rap1-peak-calls.bed -o Rap1-features.bed
#
#   For help use: findMyFeatures.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Peak calls
#   -o, --output        Features
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
require 'master_annotations'
require 'bed'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i readcount.wig -w orfs.bed -o ORFWindows.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Peak Calls (Bed)" ) { |f| options[:input] = f }
  opts.on( '-o', '--output FILE', :required, "Features (sacCer2)" ) { |f| options[:output] = f }
      
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end


puts 'Loading the MasterAnnotations DB'
loci = MasterAnnotations.to_spot_array
  
puts "\nLoading binding sites"
bs = Bed.load(options[:input])
  
puts "\nSorting MasterAnnotations entries into bound and unbound groups"
bound = SpotArray.new
loci.each do |chr,spots|
  bound[chr] ||= Array.new
  
  # Sort the spots into the two groups
  spots.each do |spot|
		# Move on if there are no binding sites on this chromosome
		next if bs[chr].nil?
	
    # If any binding sites are within this spot, put in bound group
    if bs[chr].select { |site| spot.include?(site.start) or spot.include?(site.stop) }.length > 0
      bound[chr] << spot
    end
  end
end

puts "\nWriting bound features to disk"
File.open(options[:output],'w') do |f|
	bound.each do |chr,spots|
		spots.each do |spot|
			f.puts "#{spot.id}\t{chr}\t#{spot.start}\t#{spot.stop}"
		end
	end
end