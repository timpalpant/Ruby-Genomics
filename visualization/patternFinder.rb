#!/usr/bin/env ruby1.9

# == Synopsis 
#   Find patterns in sequencing files by averaging the values in a list of Bed windows
#		Essentially matrixAligner, but keep only the average value and not the whole matrix	
#		Can also average multiple files at once to the same loci
#
# == Usage 
#   Average readcount.wig to loci specified in TSS.bed and output
#   the values to TSS-average-profile.txt
#
#   patternFinder.rb -l TSS.txt -o TSS-average-profile.txt input1.wig input2.wig
#
#   For help use: patternFinder.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -l, --loci          List of loci to align to (chromosome  start  stop  id   alignmentBase)
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
require 'wig'
require 'bed'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -l orfs.txt -o orfs-aligned.txt input1.wig input2.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-l', '--loci FILE', :required, "List of loci to align to (Bed format)" ) { |f| options[:loci] = f }
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
  	
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches? or ARGV.length < 1
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end


# Load the list of loci to align to
puts 'Loading the list of alignment loci' if ENV['DEBUG']
loci = Bed.load(options[:loci])

puts "\nComputing alignment dimensions\n" if ENV['DEBUG']
left_max = loci.collect { |chr,spots| spots.collect { |spot| (spot.value-spot.start).abs }.max }.max.to_i
right_max = loci.collect { |chr,spots| spots.collect { |spot| (spot.value-spot.stop).abs }.max }.max.to_i
# One bonus for odd/even safety
n = left_max + right_max + 1
alignment_point = left_max
puts "Average will be computed for #{n} bases from #{-alignment_point} to #{n-alignment_point}" if ENV['DEBUG']
indices = (-alignment_point..n-alignment_point).to_a

puts "\nInitializing input files\n" if ENV['DEBUG']
wigs = ARGV.map { |input_file| BigWigFile.new(input_file) }

puts "\nBeginning averaging\n" if ENV['DEBUG']
averages = Array.new
wigs.each do |wig|
	# Align and average values from all loci
	puts "\nAveraging values for: #{File.basename(wig.data_file)}" if ENV['DEBUG']
	sum = Array.new(n, 0)
	count = Array.new(n, 0)
	loci.each do |chr,spots|
		spots.each do |spot|
			begin
				values = wig.query(chr, spot.start, spot.stop)
			rescue GenomicIndexError
				next
			end
		  
		  # Locus alignment point (spot.value) should be positioned over
		  # the matrix alignment point (alignment_point)
		  n1 = alignment_point - (spot.value - spot.start).abs.to_i
		  n2 = alignment_point + (spot.value - spot.stop).abs.to_i
		  # length we are trying to insert should equal the length we are replacing
			raise "Spot is not the right length!: #{values.length} vs. #{n2-n1+1}, ({chr},#{spot})" if values.length != (n2-n1+1)
		
		  # Add values to the sum and add one to the count for those bases
      for bp in n1..n2
        sum[bp] += values[bp-n1]
        count[bp] += 1
      end
		end
	end

	puts 'Computing average' if ENV['DEBUG']
  avg = Array.new(sum.length, 'NaN')
  for i in 0...sum.length
    avg[i] = sum[i] / count[i] unless count[i] == 0
  end
	averages << avg
end

puts 'Writing averages to disk' if ENV['DEBUG']
File.open(options[:output],'w') do |f|
	for i in 0...averages[0].length
		f.puts "#{indices[i]}\t#{averages.collect { |avg| avg[i] }.join("\t")}"
	end
end
