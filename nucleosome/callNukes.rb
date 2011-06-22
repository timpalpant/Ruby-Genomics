#!/usr/bin/env ruby1.9

# == Synopsis 
#   Identifies nucleosome centers in a Wig file
#
# == Usage 
#   Identify nucleosome centers in file1.wig:
#
#   callNukes.rb -r absoluteDyads.wig -s smoothedDyads.wig -n 147 -o nucCenters.wig
#
#   For help use: callNukes.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -a, --absolute      Absolute positioning (dyads) file
#   -s, --smoothed      Smoothed positioning (dyads) file
#   -n, --nuke-size     Assumed size of a nucleosome (in base pairs)
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
require 'stats'
require 'nucleosome'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -a dyads.bw -c smoothededDyads.bw -n 147 -o output.nukes"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-a', '--absolute FILE', :required, "Absolute positioning (dyads) file" ) { |f| options[:absolute] = f }
  opts.on( '-c', '--smoothed FILE', :required, "smoothed positioning file (smoothed)" ) { |f| options[:smoothed] = f }
  options[:nuke] = 147;
  opts.on( '-n', '--nuke-size BP', "Assumed size of a nucleosome (in base pairs)" ) { |num| options[:nuke] = num.to_i }
  opts.on( '-o', '--output FILE', :required, "Output file (Nucleosome centers)" ) { |f| options[:output] = f }
      
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end


# Half nuke size for searching windows
half_nuke = options[:nuke] / 2

# Initialize WigFile files
absolute = BigWigFile.new(options[:absolute])
smoothed = BigWigFile.new(options[:smoothed])
  
# Validate that all files have the same chromosomes
absolute.chromosomes.each do |chr_id|
  raise("Files have different chromosomes!") unless smoothed.chromosomes.include?(chr_id)
end
  
File.open(options[:output],'w') do |f|
  # Header line
  f.puts NukeCalls::HEADER

  # Iterate over all chromosomes
  absolute.chromosomes.each do |chr|
    # Load the chromosome
    absolute_chr = absolute[chr]
    # Map nil values on the ends of chromosomes (from smoothing) to zero
    smoothed_chr = smoothed[chr].map { |value| value.nil? ? 0 : value }

    smoothed_chr.sort_index.reverse.each do |i|
      if smoothed_chr[i] > 0
        nuke = Nucleosome.new
        nuke.start = [0, i-half_nuke].max
        nuke.stop = [i+half_nuke, smoothed_chr.length-1].min
        nuke.dyad = i
        
        nuke.occupancy = 0
        weighted_sum, smoothed_sum = 0, 0
        for bp in nuke.start..nuke.stop
          nuke.occupancy += absolute_chr[bp]
          weighted_sum += absolute_chr[bp] * bp
          smoothed_sum += smoothed_chr[bp]
        end
        nuke.occupancy = nuke.occupancy.to_i

        if nuke.occupancy > 0
          nuke.dyad_mean = (weighted_sum / nuke.occupancy).to_i
          nuke.conditional_position = smoothed_chr[i] / smoothed_sum

          sum_of_squares = 0
          for bp in nuke.start..nuke.stop
            sum_of_squares += absolute_chr[bp] * (bp - nuke.dyad_mean)**2
          end
          nuke.dyad_stdev = Math.sqrt(sum_of_squares.to_f / nuke.occupancy)
          
          # Write nucleosome to output
          f.puts nuke.to_s
          
          # Set 147bp (nuke size) surrounding current bp on either side to 0
          # This is the region in which another nuke cannot be called
          low = Math.max(i-options[:nuke], 0)
          high = Math.min(i+options[:nuke], smoothed_chr.length-1)
          for bp in low..high
            smoothed_chr[bp] = 0
          end
        end
      end
    end
  end
end
