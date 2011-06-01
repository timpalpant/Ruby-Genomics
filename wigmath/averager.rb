#!/usr/bin/env ruby1.9

# == Synopsis 
#   Averages multiple Wig Files
#
# == Usage 
#   Average file1.wig and file2.wig:
#
#   averager.rb file1.wig file2.wig -o wig
#
#   For help use: averager.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -g, --genome        Genome assembly to use (in common/genomes/*)
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
require 'forkmanager'
require 'pickled_optparse'
require 'wig'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} file1.wig file2.wig -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  options[:step] = 200_000
  opts.on( '-c', '--step N', "Chunk size to use in base pairs (default: 200,000)" ) { |n| options[:step] = n.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
      
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
end


# Initialize the wig files to average
wigs = ARGV.map { |filename| WigFile.new(filename) }
num_files = wigs.length.to_f
# Validate their compatibility
wigs[1..-1].each do |wig|
  wigs.first.chromosomes.each do |chr|
    raise "Wig files do not have the same chromosomes" unless wig.include?(chr)
    raise "Wig file chromosome lengths are incompatible" unless wigs.first.chr_length(chr) == wig.chr_length(chr)
  end
end

# Initialize the process manager
pm = Parallel::ForkManager.new(options[:threads])

# Process each chromosome in parallel
wigs.first.chromosomes.each do |chr|
  # Run in parallel processes managed by ForkManager
  pm.start(chr) and next
  
  puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']

	# Write the chromosome fixedStep header
	File.open(options[:output]+'.'+chr, 'w') do |f|
		f.puts Wig.fixed_step(chr) + ' start=1 step=1 span=1'
	end
	
	chunk_start = 1
  chr_length = wigs.first.chr_length(chr)
	while chunk_start < chr_length
    chunk_stop = chunk_start + options[:step] - 1
    puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
    
		sum = wigs.first.query(chr, chunk_start, chunk_stop)
    wigs[1..-1].each do |wig|
      data = wig.query(chr, chunk_start, chunk_stop)
      for i in 0...data.length
        sum[i] += data[i]
      end
    end

    avg = sum.map { |value| value / num_files }

		# Write this chunk to disk
		File.open(options[:output]+'.'+chr, 'a') do |f|
			f.puts avg.map { |value| value.to_s(5) }.join("\n")
		end
		
		chunk_start = chunk_stop + 1
	end

  pm.finish(0)
end

# Wait for all of the child processes (each chromosome) to complete
pm.wait_all_children

# Iterate over the Wig file chromosome-by-chromosome
header_file = options[:output]+'.header'
File.open(header_file, 'w') do |f|
  name = "Average of #{ARGV.map { |f| File.basename(f) }.join(',')}"
  f.puts Wig.track_header(name,name)
end

# Concatenate all of the individual chromosomes into the output file
tmp_files = [header_file]
wigs.first.chromosomes.each { |chr| tmp_files << (options[:output]+'.'+chr) }
File.cat(tmp_files, options[:output])

# Delete the individual chromosome files created by each process
tmp_files.each { |filename| File.delete(filename) }