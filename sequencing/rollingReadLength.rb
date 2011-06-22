#!/usr/bin/env ruby1.9

# == Synopsis 
#   Computes the rolling read length across the genome
#
# == Usage 
#   Take BAM reads and computes the avg read length
#   over all bases in the genome
#
#   rollingReadLength.rb -i bowtie.lengths.bed -o readLengths.bw
#
#   For help use: rollingReadLength.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (Bed)
#   -g, --genome        Genome assembly to use
#   -o, --output        Output file with read length averages (BigWig)
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
require 'unix_file_utils'
require 'forkmanager'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.sam -o readlengths.bw"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file with reads (BAM)" ) { |f| options[:input] = f }
  options[:genome] = 'sacCer2'
  opts.on( '-g', '--genome ASSEMBLY', "Assembly to use (default: sacCer2)" ) { |g| options[:genome] = g }
  options[:step] = 100_000
  opts.on( '-s', '--step N', "Initial step size to use in base pairs (default: 100,000)" ) { |n| options[:step] = n.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
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


# Index the BAM file for random lookup
SAMTools.index(options[:input]) if not File.exist?(options[:input]+'.bai')

# Initialize the assembly to generate coverage on
assembly = Assembly.load(options[:genome])

# Initialize the process manager
pm = Parallel::ForkManager.new(options[:threads])


# Process each chromosome in chunks
# Each chromosome in a different parallel process
assembly.each do |chr, chr_length|
  # Run in parallel processes managed by ForkManager
  pm.start(chr) and next
  
  puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
  
  # Write the chromosome fixedStep header
  File.open(options[:output]+'.'+chr, 'w') do |f|
    f.puts Wig.fixed_step(chr) + ' start=1 step=1 span=1'
  end
  
	# Write the chromosome fixedStep header
	File.open(options[:output]+'.'+chr, 'w') do |f|
		f.puts WigFile.fixed_step(chr) + ' start=1 step=1 span=1'
	end
	
	chunk_start = 1
	while chunk_start < chr_length		
		# Allocate memory for this chunk
		chunk_stop = [chunk_start+options[:step]-1, chr_length].min
    chunk_size = chunk_stop - chunk_start + 1
    total = Array.new(chunk_size, 0)
    read_count = Array.new(chunk_size, 0)
    
    # Pad the query because SAMTools only returns reads that physically overlap the given window
    # It is possible that our 36bp reads may be physically outside the chunk window
    # but will extend into the chunk window
    query_start = [chunk_start-500, 1].max
    query_stop = [chunk_stop+500, chr_length].min
    
    # Count the number of reads for this chunk to make sure it's a reasonable number
    # Adjust the step size to an optimal size
    count = SAMTools.count(options[:input], chr, chunk_start, chunk_stop)
    puts "#{count} reads in block #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
    if count > 500_000
      options[:step] = 3*options[:step]/5
      puts "Shrinking step size - now #{options[:step]}" if ENV['DEBUG']
      redo
    elsif count < 100_000 and options[:step] < 1_000_000 and chunk_size == options[:step]
      options[:step] *= 2
      puts "Increasing step size - now #{options[:step]}" if ENV['DEBUG']
      redo
    end
    
    # Get all aligned reads for this chunk and map the dyads
    SAMTools.view(options[:input], chr, query_start, query_stop).each do |entry|      
      # Also clamp to the chunk
      low = [entry.low-chunk_start, 0].max
      high = [entry.high-chunk_start, chunk_size-1].min
      
      # Map the read coverage within the wig file
      for bp in low..high
        total[bp] += entry.length
        read_count[bp] += 1
      end
    end
    
    mean = Array.new(total.length, 0)
    for i in 0...total.length
      mean[i] = total[i].to_f / read_count[i] unless read_count[i] == 0
    end
  
    # Write this chunk to disk
    File.open(options[:output]+'.'+chr, 'a') do |f|
      f.puts mean.join("\n")
    end
    
    chunk_start = chunk_stop + 1
  end
  
  pm.finish(0)
end


# Wait for all of the child processes (each chromosome) to complete
pm.wait_all_children

# Write the Wiggle track header
header_file = options[:output]+'.header'
File.open(header_file, 'w') do |f|
	name = "Mean Read Length #{File.basename(options[:input])}"
	f.puts UCSCTools.track_header(:name => name)
end

# Concatenate all of the individual chromosomes into the output file
tmp_files = [header_file]
assembly.chromosomes.each { |chr| tmp_files << (options[:output]+'.'+chr) }
File.cat(tmp_files, options[:output])

# Delete the individual chromosome files created by each process
tmp_files.each { |filename| File.delete(filename) }

# Delete the BAM index so that it is not orphaned within Galaxy
File.delete(options[:input] + '.bai')

# Conver the output Wig file to BigWig
tmp = options[:output] + '.tmp'
Wig.to_bigwig(options[:output], tmp, assembly)
FileUtils.move(tmp, options[:output])
