#!/usr/bin/env ruby1.9

# == Synopsis 
#   Find the division of 2 Wig files
#
# == Usage 
#   Divide file2.wig into file1.wig:
#
#   divide.rb -m file1.wig -s file2.wig -o output.wig
#
#   For help use: divide.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -1, --dividend      File 1
#   -2, --divisor       File 2
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
  opts.banner = "Usage: ruby #{__FILE__} -1 file1.wig -2 file2.wig -o output.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-1', '--dividend FILE', :required, "File 1" ) { |f| options[:dividend] = f }
  opts.on( '-2', '--divisor FILE', :required, "File 2" ) { |f| options[:divisor] = f }
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

# Initialize WigFile files
dividend = WigFile.new(options[:dividend])
divisor = WigFile.new(options[:divisor])
  
# Validate that both files have the same chromosomes
puts "Validating compatibility" if ENV['DEBUG']
dividend.chromosomes.each do |chr_id|
	# TODO: Also check compatibility of start, step, span
  unless divisor.chromosomes.include?(chr_id)
		raise "Files have different/incompatible chromosomes!" 
	end
end

File.open(options[:output],'w') do |f|
  name = "#{File.basename(options[:dividend])} / #{File.basename(options[:divisor])}"
  desc = "Ratio of #{File.basename(options[:dividend])} / #{File.basename(options[:divisor])}"
  f.puts Wig.track_header(name, desc)
        
  dividend.chromosomes.each do |chr_id|
    puts "Processing chromosome #{chr_id}" if ENV['DEBUG']
    
    # Load the chromosome from both files
    dividend_chr = dividend[chr_id]
    divisor_chr = divisor[chr_id]
    
    # Check that the chromosomes have the same number of values
    raise "Chromosome #{chr_id} has different # of values! (#{dividend_chr.length} vs. #{divisor_chr.length}" if dividend_chr.length != divisor_chr.length
    
    # Compute the ratio for all values in the chromosome
		ratio = Chromosome.new(dividend_chr.length, dividend_chr.start, dividend_chr.step, dividend_chr.span)
    ratio[0..-1] = (dividend_chr / divisor_chr)
    
    # Write to file
    f.puts Wig.fixed_step(chr_id, ratio)
    f.puts ratio
  end
end