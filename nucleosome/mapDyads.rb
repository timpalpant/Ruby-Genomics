#!/usr/bin/env ruby1.9

# == Synopsis 
#   Maps the density of read centers
#
# == Usage 
#   Take paired-end SAM reads and count the number of read centers
#   at each loci. (Alternative to baseAlignCounts)
#
#   mapDyads.rb -i bowtie.bed -o dyads.wig
#
#   For help use: mapDyads.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (BAM)
#		-l, --length				Mononucleosome length (for offset)
#   -o, --output        Output file with read center density (Wig)
#   -g, --genome        Genome assembly to use (in common/assemblies/*)
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
require 'wig'
require 'sam'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.bam -o dyads.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input file with reads (BAM)" ) { |f| options[:input] = f }
  opts.on( '-l', '--length N', "Mononucleosome length (default: read length)" ) { |n| options[:length] = n.to_i }
	options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome NAME', "Genome assembly (default: sacCer2)" ) { |name| options[:genome] = name }
	options[:step] = 1_000_000
  opts.on( '-s', '--step N', "Step size to use in base pairs (default: 500,000)" ) { |n| options[:step] = n.to_i }
	opts.on( '-o', '--output FILE', :required, "Output file (Wig)" ) { |f| options[:output] = f }
      
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end	
end

# Warning if using manual offset
if options[:length]
	offset = options[:length] / 2
	puts "Using fixed offset of #{offset} from read starts (5' end)" if ENV['DEBUG']
end

# Load the genome assembly
assembly = Assembly.load(options[:genome])

# Initialize the BAM file
bam = BAMFile.new(options[:input])

# Write the Wiggle track header
File.open(options[:output], 'w') do |f|
	name = "Mapped starts #{File.basename(options[:input])}"
	f.puts Wig.track_header(name, name)
end

unmapped = 0
# Process each chromosome in options[:step] bp chunks
assembly.each do |chr, chr_length|
	print "\nProcessing chromosome #{chr}" if ENV['DEBUG']

	# Write the chromosome fixedStep header
	File.open(options[:output], 'a') do |f|
		f.puts Wig.fixed_step(chr) + ' start=1 step=1 span=1'
	end
	
	chunk_start = 0
	while chunk_start < chr_length
		# Some indication of progress
		print '.' if ENV['DEBUG']
		
		# Allocate memory for this chunk
		chunk_size = [options[:step], chr_length-chunk_start].min
		mapped_starts = Array.new(chunk_size, 0)
	
		# Get all aligned reads for this chunk and map the dyads
		bam.foreach(chr, chunk_start, chunk_start + options[:step] - 1) do |read|
			begin
				center = if options[:length]
					if read.watson?
						read.start + offset
					else
						read.start - offset
					end
				else
					read.center
				end
				
				mapped_starts[center-start] += 1
			rescue
				unmapped += 1
			end
		end
		
		# Write this chunk to disk
		File.open(options[:output], 'a') do |f|
			f.puts mapped_starts.join("\n")
		end
		
		chunk_start += options[:step]
	end
	
	# Force GC after every chromosome (bio-samtools is leaky?)
	10.times { GC.start }
end

puts "WARN: #{unmapped} unmapped reads" if unmapped > 0

# Close the BAM file
bam.close