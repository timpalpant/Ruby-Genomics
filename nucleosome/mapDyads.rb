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
require 'samtools'

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
	options[:step] = 100_000
  opts.on( '-s', '--step N', "Initial step size to use in base pairs (default: 100,000)" ) { |n| options[:step] = n.to_i }
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

# Index the BAM file for random lookup
SAMTools.index(options[:input])

# Load the genome assembly
assembly = Assembly.load(options[:genome])

# Write the Wiggle track header
File.open(options[:output], 'w') do |f|
	name = "Mapped starts #{File.basename(options[:input])}"
	f.puts Wig.track_header(name, name)
end

unmapped = 0
# Process each chromosome in options[:step] bp chunks
assembly.each do |chr, chr_length|
	puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']

	# Write the chromosome fixedStep header
	File.open(options[:output], 'a') do |f|
		f.puts Wig.fixed_step(chr) + ' start=1 step=1 span=1'
	end
	
	chunk_start = 0
	while chunk_start < chr_length		
		# Allocate memory for this chunk
		chunk_size = [options[:step], chr_length-chunk_start].min
		mapped_starts = Vector::Int[chunk_size] #Array.new(chunk_size, 0)
    chunk_stop = chunk_start + chunk_size - 1
    
    # Count the number of reads for this chunk to make sure it's a reasonable number
    # Adjust the step size to an optimal size
    count = SAMTools.count(options[:input], chr, chunk_start, chunk_stop)
    puts "#{count} reads in block #{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
    if count > 1_000_000
      options[:step] =4*options[:step]/5
      puts "Shrinking step size - now #{options[:step]}" if ENV['DEBUG']
      redo
    elsif count < 100_000 and options[:step] < 5_000_000 and chunk_size == options[:step]
      options[:step] *= 2
      puts "Increasing step size - now #{options[:step]}" if ENV['DEBUG']
      redo
    end
	
		# Get all aligned reads for this chunk and map the dyads
		SAMTools.view(options[:input], chr, chunk_start, chunk_stop).each do |read|
			center = if options[:length]
				if read.watson?
					read.start + offset
				else
					read.start - offset
				end
			else
				read.center
			end
				
			begin
				mapped_starts[center-chunk_start] += 1 if chunk_start <= center and center < chunk_start+chunk_size
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
end

puts "WARN: #{unmapped} unmapped reads" if unmapped > 0

# Delete the BAM index so that it is not orphaned within Galaxy
File.delete(options[:input] + '.bai')