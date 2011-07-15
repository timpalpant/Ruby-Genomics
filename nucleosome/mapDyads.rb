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
#   -l, --length        Mononucleosome length (for offset)
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
require 'utils/parallelizer'
require 'bio-genomic-file'
require 'utils/unix'
require 'reference_assembly'
include Bio

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
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
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
assembly = ReferenceAssembly.load(options[:genome])

# Process each chromosome in chunks
# Each chromosome in a different parallel process
unmapped_counts = nil
BAMFile.open(options[:input]) do |bam|
  unmapped_counts = assembly.p_map(:in_processes => options[:threads]) do |chr, chr_length|
    puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
    unmapped = 0

    # Write the chromosome fixedStep header
    File.open(options[:output]+'.'+chr, 'w') do |f|
      f.puts "fixedStep chrom=#{chr} start=1 step=1"
    end
    
    chunk_start = 1
    while chunk_start < chr_length		
      # Allocate memory for this chunk
      chunk_stop = [chunk_start+options[:step]-1, chr_length].min
      chunk_size = chunk_stop - chunk_start + 1
      mapped_starts = Array.new(chunk_size, 0)
      
      # Pad the query because SAMTools only returns reads that physically overlap the given window
      # It is likely that our 36bp reads may be physically outside the chunk window
      # but the dyad will be within the chunk window (since the DNA strands are actually ~150bp)
      query_start = [chunk_start-500, 1].max
      query_stop = [chunk_stop+500, chr_length].min
      
      # Get all aligned reads for this chunk and map the dyads
      bam.each_read(chr, query_start, query_stop) do |read|
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
          mapped_starts[center-chunk_start] += 1 if chunk_start <= center and center <= chunk_stop
        rescue
          unmapped += 1
        end
      end
      
      # Write this chunk to disk
      File.open(options[:output]+'.'+chr, 'a') do |f|
        f.puts mapped_starts.join("\n")
      end
      
      chunk_start = chunk_stop + 1
    end
    
    # Send the number of unmapped reads on this chromosome back to the parent process
    puts "#{unmapped} unmapped dyads on chromosome #{chr}" if unmapped > 0 and ENV['DEBUG']
    unmapped
  end
end

total_unmapped = unmapped_counts.sum
puts "WARN: #{total_unmapped} unmapped dyads" if total_unmapped > 0

# Write the Wiggle track header
header_file = options[:output]+'.header'
File.open(header_file, 'w') do |f|
  name = "Mapped Dyads #{File.basename(options[:input])}"
  f.puts Utils::UCSC::TrackHeader.new(:name => name)
end

# Concatenate all of the individual chromosomes into the output file
tmp_files = [header_file]
assembly.chromosomes.each { |chr| tmp_files << (options[:output]+'.'+chr) }
File.cat(tmp_files, options[:output])

# Delete the individual chromosome files created by each process
tmp_files.each { |filename| File.delete(filename) }

# Convert the output to BigWig
tmp = options[:output] + '.tmp'
TextWigFile.to_bigwig(options[:output], tmp, assembly)
FileUtils.move(tmp, options[:output])
