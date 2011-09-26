#!/usr/bin/env ruby1.9

# == Synopsis 
#   Aligns sequencing values for a heatmap matrix
#
# == Usage 
#   Align readcount.wig to loci specified in orfs.txt and output
#   the values to orfs-aligned.txt
#
#   matrixAligner.rb -i readcounts.bw -l orfs.bed -o orfs-aligned.matrix.txt
#
#   For help use: matrixAligner.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file to align values (BigWig)
#   -l, --loci          List of loci to align to (chromosome  start  stop  id   alignmentBase)
#   -o, --output        Output file (matrix with dimensions #(loci) x max(stop-start)
#   -m, --max           Maximum allowed length for a row (entries outside this value will be truncated)
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
require 'bio-genomic-file'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i readcount.bw -l orfs.txt -o orfs-aligned.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file to align values (BigWig)" ) { |f| options[:input] = f }
  opts.on( '-l', '--loci FILE', :required, "List of loci to align to (Bed format)" ) { |f| options[:loci] = f }
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
  opts.on( '-m', '--max N', "Maximum allowed row length" ) { |n| options[:max] = n.to_i }
  opts.on( '-s', '--markers N', "Draw markers every N base pairs" ) { |n| options[:ladder] = n.to_i }
      
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Set to whatever empty matrix values should be (e.g. NaN, -, '', etc.)
NA_PLACEHOLDER = '-'

# Load the list of loci to align to
loci = BedFile.open(options[:loci]).to_a

# Validation and default alignment points
loci.each do |spot|
  spot.start = 1 if spot.start < 1
  spot.stop = 1 if spot.stop < 1
  spot.value = spot.start if spot.value.nil? or spot.value < spot.low or spot.value > spot.high
end

m = loci.length
left_max = loci.collect { |spot| (spot.value-spot.start).abs }.max.to_i
right_max = loci.collect { |spot| (spot.value-spot.stop).abs }.max.to_i
# One bonus for odd/even safety
n = left_max + right_max + 1
alignment_point = left_max
puts "Intervals aligned into: #{m}x#{n} matrix"
puts "Alignment point: #{alignment_point}"
left_bound = 0
right_bound = n-1
if not options[:max].nil? and options[:max] < n
  puts "Truncated to #{m}x#{options[:max]}"
  left_align_distance = alignment_point
  right_align_distance = n - alignment_point - 1
  half_max = options[:max] / 2
  
  # If there enough room to center the alignment point within the cutoff row length, do that
  if half_max < left_align_distance and half_max < right_align_distance
    left_bound = alignment_point - half_max
    right_bound = alignment_point + half_max
  else
    if left_align_distance <= right_align_distance
      right_bound = options[:max]
    else
      left_bound = n - options[:max]
    end
  end
end

# Construct markers for the top and bottom
if options.include?(:ladder)
  marker_line = Array.new(right_bound-left_bound+1, NA_PLACEHOLDER)
  marker = alignment_point % options[:ladder]
  while marker < right_bound-left_bound+1
    marker_line[marker] = '1e30'
    marker += options[:ladder]
  end
  marker_line ="Marker\t" +  marker_line.join("\t")
end


# Initialize the Wig file
puts "Initializing BigWig file" if ENV['DEBUG']
wig = WigFile.autodetect(options[:input])

# Align values for each locus around the alignment point
skipped = 0
puts "Aligning matrix and writing to disk" if ENV['DEBUG']
File.open(options[:output],'w') do |f|
  # Write a header (required by matrix2png)
  f.puts "ID\t" + (left_bound-alignment_point..right_bound-alignment_point).to_a.join("\t")
  # Add markers
  10.times { f.puts marker_line } if options.include?(:ladder)
  
  # Write the data for each row
  loci.each_with_index do |spot,i|
    # Get the data for this interval from the wig file
    begin
      result = wig.query(spot.chr, spot.low, spot.high)
      values = (spot.low..spot.high).map { |bp| result.get(bp) }
      values.reverse! if spot.crick?
    rescue
      puts "Skipping spot (#{spot.id},#{spot.chr}:#{spot.start}-#{spot.stop},#{spot.value}) because data could not be retrieved" if ENV['DEBUG']
      skipped += 1
      next
    end
      
    # Locus alignment point (spot.value) should be positioned over
    # the matrix alignment point (alignment_point)
    n1 = alignment_point - (spot.value-spot.start).abs.to_i
    n2 = alignment_point + (spot.value-spot.stop).abs.to_i
    # length we are trying to insert should equal the length we are replacing
    raise "Spot is not the right length!: #{values.length} vs. #{n2-n1+1}, (#{spot})" if values.length != (n2-n1+1)
      
    entry = Array.new(n, NA_PLACEHOLDER)
    entry[n1..n2] = values.map { |v| v ? v : NA_PLACEHOLDER }
    # Total length should be the matrix width to avoid irregular matrices
    raise "Entry is not the right length!: #{entry.length} vs. #{n}, ({chr},#{spot})" if entry.length != n
      
    # Save to the output file
    id = spot.id ? spot.id : "Row #{i+1}"
    f.puts id + "\t" + entry[left_bound..right_bound].join("\t")
  end
  
  # Add markers at the bottom
  10.times { f.puts marker_line } if options.include?(:ladder)
end

puts "#{skipped} spots skipped"
