#!/usr/bin/env ruby1.9

# == Synopsis 
#   Maps the density of reads
#
# == Usage 
#   Take BAM reads and count the number of reads overlapping
#   at each loci (occupancy).
#
#   baseAlignCounts.rb -i bowtie.bam -o occupancy.bw
#
#   For help use: baseAlignCounts.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (BAM)
#   -o, --output        Output file with occupancy (BigWig)
#   -g, --genome        Genome assembly to use
#   -x, --extend        In silico artificial extension
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
require 'utils/parallelizer'
require 'pickled_optparse'
require 'bio-genomic-file'
require 'utils/unix'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.bam -o occupancy.bw"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-i', '--input FILE', :required, "Input file with reads (BAM)" ) { |f| options[:input] = f }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome NAME', "Genome assembly (default sacCer2)" ) { |name| options[:genome] = name }
  options[:step] = 100_000
  opts.on( '-s', '--step N', "Initial step size to use in base pairs (default: 100,000)" ) { |n| options[:step] = n.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-x', '--extend N', "In silico extension (default: read length)") { |n| options[:x] = n.to_i }
  opts.on( '-o', '--output FILE', :required, "Output file (BigWig)" ) { |f| options[:output] = f }
      
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Pad queries to SAMTools
if options[:x] and options[:x] > 0
  padding = options[:x] + 1
else
  padding = 500
end

# Initialize the assembly to generate coverage on
assembly = Assembly.load(options[:genome])

# Process each chromosome in chunks
# Each chromosome in a different parallel process
unmapped_counts = nil
BAMFile.open(options[:input]) do |bam|
  unmapped_counts = assembly.p_map(:in_processes => options[:threads]) do |chr, chr_length|
    puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
    unmapped = 0

    # Write the chromosome fixedStep header
    File.open(options[:output]+'.'+chr, 'w') do |f|
      f.puts Contig.new(0, chr, 1, 1, 1).to_s
    end
    
    chunk_start = 1
    while chunk_start < chr_length		
      # Allocate memory for this chunk
      chunk_stop = [chunk_start+options[:step]-1, chr_length].min
      chunk_size = chunk_stop - chunk_start + 1
      occupancy = Array.new(chunk_size, 0)
      
      # Pad the query because SAMTools only returns reads that physically overlap the given window
      # It is possible that our 36bp reads may be physically outside the chunk window
      # but will extend into the chunk window
      query_start = [chunk_start-padding, 1].max
      query_stop = [chunk_stop+padding, chr_length].min
    
      # Get all aligned reads for this chunk and map the dyads
      bam.each(chr, query_start, query_stop) do |entry|
        # Calculate the read stop based on specified in silico extension
        # or the read length
        stop = if not options[:x].nil? and options[:x] != 0
          strand = (entry.watson? ? 1 : -1)
          entry.start + strand*options[:x] - 1
        else
          entry.stop
        end

        # Get the high and low read coordinates, and clamp to the ends of the chromosome
        low = [1, [entry.start, stop].min].max
        high = [[entry.start, stop].max, chr_length].min
        
        # Also clamp to the chunk
        low = [low-chunk_start, 0].max
        high = [high-chunk_start, chunk_size-1].min
      
        # Map the read coverage within the wig file
        begin
          for bp in low..high
            occupancy[bp] += 1
          end
        rescue
          unmapped += 1
        end
      end
      
      # Write this chunk to disk
      File.open(options[:output]+'.'+chr, 'a') do |f|
        f.puts occupancy.join("\n")
      end
      
      chunk_start = chunk_stop + 1
    end
    
    # Send the number of unmapped reads on this chromosome back to the parent process
    puts "#{unmapped} unmapped reads on chromosome #{chr}" if unmapped > 0 and ENV['DEBUG']
    unmapped
  end
end

total_unmapped = unmapped_counts.sum
puts "WARN: #{total_unmapped} unmapped reads" if total_unmapped > 0

# Write the Wiggle track header
header_file = options[:output]+'.header'
File.open(header_file, 'w') do |f|
  name = "Mapped Coverage #{File.basename(options[:input])}"
  f.puts UCSCTrackHeader.new(:name => name)
end

# Concatenate all of the individual chromosomes into the output file
tmp_files = [header_file]
assembly.chromosomes.each { |chr| tmp_files << (options[:output]+'.'+chr) }
File.cat(tmp_files, options[:output])

# Delete the individual chromosome files created by each process
tmp_files.each { |filename| File.delete(filename) }

# Conver the output Wig file to BigWig
tmp = options[:output] + '.tmp'
Wig.to_bigwig(options[:output], tmp, assembly)
FileUtils.move(tmp, options[:output])
