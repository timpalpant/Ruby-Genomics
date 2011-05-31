#!/usr/bin/env ruby1.9

# == Synopsis 
#   Gaussian smooths Wig Files
#
# == Usage 
#   Smooth file1.wig:
#
#   smoother.rb -i file1.wig -o file1.smoothed.wig
#
#   For help use: smoother.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file
#   -s, --stdev         Standard deviation of the Gaussian
#   -w, --window        Window size
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
require 'fixed_precision'
require 'wig'
require 'stats'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i file1.wig -o file1.smoothed.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file" ) { |f| options[:input] = f }
  options[:sdev] = 10
  opts.on( '-s', '--sdev NUM', "Standard deviation of the Gaussian in base pairs (default 10)" ) { |num| options[:sdev] = num.to_i }
  options[:window_size] = 3
  opts.on( '-w', '--window NUM', "Number of standard deviations +/- to make a window (default 3)" ) { |num| options[:window_size] = num.to_i }
  options[:step] = 100_000
  opts.on( '-s', '--step N', "Chunk size to use in base pairs (default: 100,000)" ) { |n| options[:step] = n.to_i }
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

# Gaussian smoothing requires padding of half_window on either end
padding = options[:sdev] * options[:window_size] / 2
  
# Initialize the wig files to smooth
wig = WigFile.new(options[:input])

# Initialize the process manager
pm = Parallel::ForkManager.new(options[:threads])


# Process each chromosome in chunks
# Each chromosome in a different parallel process
wig.chromosomes.each do |chr|
  # Run in parallel processes managed by ForkManager
  pm.start(chr) and next
  
  puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']

	# Write the chromosome fixedStep header
	File.open(options[:output]+'.'+chr, 'w') do |f|
		f.puts Wig.fixed_step(chr) + ' start=1 step=1 span=1'
	end
	
	chunk_start = 1
  chr_length = wig.chr_length(chr)
	while chunk_start < chr_length
    chunk_stop = chunk_start + options[:step] - 1
    
		chunk = wig.query(chr, chunk_start-padding, chunk_stop+padding)
		smoothed = chunk.gaussian_smooth(options[:sdev], options[:window_size])[padding..-padding]
    
		# Write this chunk to disk
		File.open(options[:output]+'.'+chr, 'a') do |f|
			f.puts smoothed.map { |value| value.to_s(5) }.join("\n")
		end
		
		chunk_start = chunk_stop + 1
	end
end


# Wait for all of the child processes (each chromosome) to complete
pm.wait_all_children

# Iterate over the Wig file chromosome-by-chromosome
header_file = options[:output]+'.header'
File.open(header_file, 'w') do |f|
  name = "Smoothed #{File.basename(options[:input])}"
  f.puts Wig.track_header(name,name)
end

# Concatenate all of the individual chromosomes into the output file
tmp_files = [header_file]
wig.chromosomes.each { |chr| tmp_files << (options[:output]+'.'+chr) }
File.cat(tmp_files, options[:output])

# Delete the individual chromosome files created by each process
tmp_files.each { |filename| File.delete(filename) }