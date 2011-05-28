#!/usr/bin/env ruby1.9

# == Synopsis 
#   Aligns sequencing values for a heatmap matrix
#
# == Usage 
#   Generate 2 heatmaps of the windows in Intergenics.bed
#   One for bound and one for unbound using the binding sites in Nhp6a-bs.bed
#
#   boundVsUnbound.rb -i readcount.wig -b Nhp6a-bs.bed -l Intergenics.bed -o Nhp6A-Intergenics.txt
#
#   For help use: boundVsUnbound.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with heatmap data (Wig)
#   -l, --loci          List of windows to align (chr  start  stop  id   alignmentBase)
#   -b, --binding       Binding sites (chr locus)
#   -o, --output        Output files (matrix with dimensions #(loci) x max(stop-start)
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
require 'wig'
require 'bed'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i readcount.wig -l orfs.txt -o orfs-aligned.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file to align values (Wig)" ) { |f| options[:input] = f }
  opts.on( '-l', '--loci FILE', :required, "List of loci to align to (Bed format)" ) { |f| options[:loci] = f }
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
	opts.on( '-b', '--binding FILE', :required, "Binding sites (chr locus)" ) { |f| options[:bs] = f }
      
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
  
  # Generate output names for bound and unbound loci
  filename = options[:output].split('.')
  basename = filename[0..-2].join('.')
  extension = filename.last
  options[:bound] = basename + '.bound.' + extension
  options[:unbound] = basename + '.unbound.' + extension
end


puts 'Loading the list of alignment loci'
loci = Bed.load(options[:loci])
  
puts "\nLoading binding sites"
bs = Bed.load(options[:bs])
  
puts "\nSorting loci into bound and unbound groups"
bound = SpotArray.new
unbound = SpotArray.new
loci.each do |chr,spots|
  bound[chr] ||= Array.new
  unbound[chr] ||= Array.new
  
  # Sort the spots into the two groups
  spots.each do |spot|
    # If any binding sites are within this spot, put in bound group
    if bs[chr].select { |site| spot.include?(site.start) }.length > 0
      bound[chr] << spot
    # Otherwise, put in unbound group
    else
      unbound[chr] << spot
    end
  end
end

puts "\nComputing matrix dimensions"
m1 = bound.num_spots
m2 = unbound.num_spots
left_max = loci.collect { |chr,spots| spots.collect { |spot| (spot.value-spot.start).abs }.max }.max.to_i
right_max = loci.collect { |chr,spots| spots.collect { |spot| (spot.value-spot.stop).abs }.max }.max.to_i
# One bonus for odd/even safety
n = left_max + right_max + 1
alignment_point = left_max
puts "Will create a #{m1}x#{n} matrix of bound loci"
puts "Will create a #{m2}x#{n} matrix of unbound loci"
puts "Alignment point is #{alignment_point}"
  
# Load the Wig file with data
wig = Wig.load(options[:input])

# Generate a heatmap
def heatmap(filename, loci, data, n, alignment_point)
  # Generate a heatmap for both bound windows
  puts "\nGenerating heatmap #{File.basename(filename)}"
  puts 'Aligning values in the matrix'
  row_num = 1
  File.open(filename,'w') do |f|
    File.open(filename+'.key','w') do |g|
      loci.each do |chr,spots|
        puts "Processing chromosome #{chr}"
        spots.select { |spot| data[chr].include?(spot.low,spot.high) }.each do |spot|
          values = data[chr][spot.low..spot.high]
          # reverse if encoded on the minus/Crick strand (right-to-left)
          values.reverse! if spot.crick?
          
          # Locus alignment point (spot.value) should be positioned over
          # the matrix alignment point (alignment_point)
          n1 = alignment_point - (spot.value - spot.start).abs.to_i
          n2 = alignment_point + (spot.value - spot.stop).abs.to_i
          # Length we are trying to insert should equal the length we are replacing
          raise "Spot is not the right length!: #{values.length} vs. #{n2-n1+1}, (#{chr},#{spot})" if values.length != (n2-n1+1)
          
          # Make a new row in the heatmap matrix and insert our values
          entry = Array.new(n, 'NaN')
          entry[n1..n2] = values.to_a
          # Total length should be the matrix width to avoid irregular matrices
          raise "Entry is not the right length!: #{entry.length} vs. #{n}, (#{chr},#{spot})" if entry.length != n
          f.puts entry.join("\t")
          g.puts "#{row_num}\t#{spot.id}\t#{chr}\t#{spot.start}\t#{spot.stop}"
          row_num += 1
        end
      end
    end
  end
  puts 'complete!'
end
  
# Generate bound and unbound heatmaps
heatmap(options[:bound], bound, wig, n, alignment_point)
heatmap(options[:unbound], unbound, wig, n, alignment_point)