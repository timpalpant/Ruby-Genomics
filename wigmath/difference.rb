#!/usr/bin/env ruby1.9

# == Synopsis 
#   Find the difference between 2 Wig files
#
# == Usage 
#   Subtract file2.wig from file1.wig:
#
#   difference.rb -m file1.wig -s file2.wig -o output.wig
#
#   For help use: difference.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -m, --minuend       File 1
#   -s, --subtrahend    File 2
#   -p, --percent       Output the percent difference
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

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -m file1.wig -s file2.wig -p -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  opts.on( '-m', '--minuend FILE', :required, "File 1" ) { |f| options[:minuend] = f }
  opts.on( '-s', '--subtrahend FILE', :required, "File 2" ) { |f| options[:subtrahend] = f }
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

global_start = Time.now

# Initialize WigFile files
minuend = WigFile.new(options[:minuend])
subtrahend = WigFile.new(options[:subtrahend])
  
# Validate that both files have the same chromosomes
puts "Validating compatibility"
minuend.chromosomes.each do |chr_id|
  raise "Files have different chromosomes!" unless subtrahend.chromosomes.include?(chr_id)
end

File.open(options[:output],'w') do |f|
  name = "#{File.basename(options[:minuend])} - #{File.basename(options[:subtrahend])}"
  desc = "Difference of #{File.basename(options[:minuend])} - #{File.basename(options[:subtrahend])}"
  f.puts Wig.track_header(name, desc)
        
  minuend.chromosomes.each do |chr_id|
    puts "Processing chromosome #{chr_id}"
    times = []
    
    # Load the chromosome from both files
    start = Time.now
    minuend_chr = minuend[chr_id]
    subtrahend_chr = subtrahend[chr_id]
    stop = Time.now
    times << (stop-start)
    
    # Check that the chromosomes have the same number of values
    raise("Chromosome #{chr_id} has different # of values! (#{minuend_chr.length} vs. #{subtrahend_chr.length}") if minuend_chr.length != subtrahend_chr.length
    
    # Compute the difference for all values in the chromosome
    start = Time.now
    difference = Chromosome.new(minuend_chr.length)
    difference.data = minuend_chr - subtrahend_chr
    stop = Time.now
    times << (stop-start)
    
    # Write to file
    start = Time.now
    f.puts Wig.fixed_step(chr_id, difference)
    f.puts difference
    stop = Time.now
    times << (stop-start)
    
    puts times.join("\t")
  end
end

global_stop = Time.now
puts "Time elapsed: #{global_stop-global_start}"