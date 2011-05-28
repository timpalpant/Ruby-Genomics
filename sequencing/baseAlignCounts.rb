#!/usr/bin/env ruby1.9

# == Synopsis 
#   Maps the density of reads
#
# == Usage 
#   Take SAM reads and count the number of reads overlapping
#   at each loci (occupancy).
#
#   baseAlignCounts.rb -i bowtie.sam -o occupancy.wig
#
#   For help use: baseAlignCounts.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (SAM)
#   -o, --output        Output file with occupancy (Wig)
#   -g, --genome        Genome assembly to use
#		-x, --extend				In silico artificial extension
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
require 'sam'
require 'math_utils'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.sam -o occupancy.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input file with reads (SAM)" ) { |f| options[:input] = f }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome NAME', "Genome assembly (default sacCer2)" ) { |name| options[:genome] = name }
  opts.on( '-o', '--output FILE', :required, "Output file (Wig)" ) { |f| options[:output] = f }
	opts.on( '-x', '--extend N', "In silico extension (default: read length)") { |n| options[:x] = n.to_i }
      
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end

# Initialize a new wig file to hold the coverage data
a = Assembly.load(options[:genome])
wig = SingleBPData.for_assembly(a)

unmapped = 0
# Iterate over the reads in the SAM file and map their coverage
SAMFile.foreach_read(options[:input]) do |entry|
	begin
		# Calculate the read stop based on specified in silico extension
		# or the read length
		stop = if not options[:x].nil? and options[:x] != 0
			strand = (entry.watson? ? 1 : -1)
			entry.start + strand*options[:x] - 1
		else
			entry.stop
		end

		# Get the high and low read coordinates, and clamp to the ends of the chromosome
		low = Math.max(1, Math.min(entry.start, stop))
		high = Math.min(Math.max(entry.start, stop), wig[entry.chr].length)
	
		# Map the read coverage within the wig file
		wig[entry.chr][low-1..high-1] += 1
	rescue
		unmapped += 1
	end
end

# Warn if there were reads that broke for some reason
# (off the end of a chromosome, 2-micron, etc.)
puts "WARN: #{unmapped} unmapped reads" if unmapped > 0

# Write the newly created Wig file to disk
wig.to_wig(options[:output])