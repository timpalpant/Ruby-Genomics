#!/usr/bin/env ruby1.9

# == Synopsis 
#   Randomly select sequencing reads
#
# == Usage 
#   Randomly select 100 reads in bowtie.bed
#
#   readSelector.rb -i bowtie.bed -n 100 -o bowtie.subset.bed
#
#   For help use: readSelector.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with mapped reads (Bed)
#   -n, --number        Number of reads to randomly select
#   -o, --output        Output file with subset of mapped reads (Bed)
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
require 'utils/unix'


# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i reads.bed -o selectedreads.bed"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file with reads (Bed)" ) { |f| options[:input] = f }
  opts.on( '-n', '--number N', :required, "Number of reads to randomly select" ) { |n| options[:num] = n.to_i }
  opts.on( '-o', '--output FILE', :required, "Output file with dumped read lengths (txt)" ) { |f| options[:output] = f }
      
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end


num_lines = File.num_lines(options[:input])
raise "Number of lines in file is less than number to be selected! Exiting..." if num_lines <= options[:num]

# If we are selecting more than half of the lines, invert our question
if options[:num] > num_lines / 2
  rejecting = true
  n = num_lines - options[:num]
else
  selecting = true
  n = options[:num]
end
    
selected = (0...num_lines).to_a.sample(n).sort
  
output_line = current_reject = 0
File.open(options[:output],'w') do |output|
  line_num = 0
  File.foreach(options[:input]) do |line|
    if output_line < options[:num] and current_reject < n
      if (selecting and selected[output_line] == line_num) or (rejecting and selected[current_reject] != line_num)
        output.write line
        output_line += 1
        puts output_line if output_line % 1_000_000 == 0 and ENV['DEBUG']
      else
        current_reject += 1
      end
    end
    
    line_num += 1
  end
end
