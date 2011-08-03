#!/usr/bin/env ruby1.9

# == Synopsis 
#   Identifies nucleosome centers in a Wig file
#
# == Usage 
#   Identify nucleosome centers in file1.wig:
#
#   callNukes.rb -r dyadsDyads.wig -s smoothedDyads.wig -n 147 -o nucCenters.wig
#
#   For help use: callNukes.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -a, --dyads      dyads positioning (dyads) file
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
require 'bio-genomic-file'
require 'stats'
require 'reference_assembly'
require 'narray'
include Bio

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
  opts.on( '-d', '--dyads FILE', :required, "dyads positioning (dyads) file" ) { |f| options[:dyads] = f }
  opts.on( '-s', '--smoothed FILE', :required, "smoothed positioning file (smoothed)" ) { |f| options[:smoothed] = f }
  opts.on( '-g', '--genome ASSEMBLY', :required, "Genome assembly" ) { |g| options[:genome] = g }
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
CHUNK_SIZE = 500_000

# Initialize the reference assembly
assembly = ReferenceAssembly.load(options[:genome])

# Initialize the Wig files
dyads = WigFile.autodetect(options[:dyads])
smoothed = WigFile.autodetect(options[:smoothed])

puts "Beginning nucleosome-calling" if ENV['DEBUG']
File.open(options[:output], 'w') do |f|
  # Header line
  f.puts NukeCallsFile::HEADER

  # Iterate over all chromosomes
  dyads.chromosomes.each do |chr|
    puts "Processing chromosome #{chr}" if ENV['DEBUG']
    chr_length = assembly[chr]
    
    # Load the entire chromosome as NArrays (for memory efficiency)
    puts "...loading data" if ENV['DEBUG']
    dyads_chr = NArray.sint(chr_length)
    smoothed_chr = NArray.sfloat(chr_length)
    
    start = 1
    while start < chr_length
      stop = [start+CHUNK_SIZE, chr_length].min
      
      dyads.query(chr, start, stop).each do |bp, value|
        dyads_chr[bp-1] = value
      end
      
      smoothed.query(chr, start, stop).each do |bp, value|
        smoothed_chr[bp-1] = value
      end
      
      start = stop + 1
    end

    puts "...sorting" if ENV['DEBUG']
    sorted_indices = smoothed_chr.sort_index
    puts "...calling nucleosomes" if ENV['DEBUG']
    for j in 1..chr_length
      # Iterate in descending order
      i = sorted_indices[chr_length-j]
      if smoothed_chr[i] > 0
        nuke = Genomics::Nucleosome.new
        nuke.chr = chr
        nuke.start = [0, i-half_nuke].max
        nuke.stop = [i+half_nuke, smoothed_chr.length-1].min
        nuke.dyad = i
        
        nuke.occupancy = 0
        weighted_sum, smoothed_sum, sum_of_squares = 0, 0, 0
        for bp in nuke.start..nuke.stop
          nuke.occupancy += dyads_chr[bp]
          weighted_sum += dyads_chr[bp] * bp
          smoothed_sum += smoothed_chr[bp]
          sum_of_squares += dyads_chr[bp] * bp**2
        end
        nuke.occupancy = nuke.occupancy.to_i

        # If we found a nucleosome, compute its dyad / fuzziness
        if nuke.occupancy > 0
          nuke.dyad_mean = (weighted_sum / nuke.occupancy).to_i
          nuke.conditional_position = smoothed_chr[i] / smoothed_sum
          variance = (sum_of_squares - weighted_sum*nuke.dyad_mean) / nuke.occupancy
          nuke.dyad_stdev = Math.sqrt(variance)
          
          # Write nucleosome to output
          f.puts nuke.to_s
          
          # Set 147bp (nuke size) surrounding current bp on either side to 0
          # This is the region in which another nuke cannot be called
          low = [i-options[:nuke], 0].max
          high = [i+options[:nuke], smoothed_chr.length-1].min
          for bp in low..high
            smoothed_chr[bp] = 0
          end
        end
      end
    end
  end
end
