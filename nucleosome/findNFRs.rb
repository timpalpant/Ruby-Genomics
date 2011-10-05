#!/usr/bin/env ruby1.9

# == Synopsis 
#   Finds the NFR by locating the +1 and -1 nucleosomes around a coordinate
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Nuke calls
#   -l, --loci          Loci to search for nukes around
#   -o, --output        NFR coordinates
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
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i nukes.txt -o nuke-distances.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file with nuke calls" ) { |f| options[:input] = f }
  opts.on( '-l', '--loci FILE', :required, "Windows file to find +1 nucleosomes (Bed format)" ) { |f| options[:loci] = f }
  opts.on( '-o', '--output FILE', :required, "Output file (nuke calls for loci)" ) { |f| options[:output] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end


skipped, invalid_coordinates = 0, 0
File.open(options[:output], 'w') do |f|
  NukeCallsFile.open(options[:input]) do |calls|    
    BedFile.foreach(options[:loci]) do |spot|
      if not calls.chromosomes.include?(spot.chr)
        skipped += 1
        next
      end
      
      # Get the +1 and -1 nucleosomes upstream and downstream of the TSS (within +/- 1kb)
      p1, m1 = nil, nil
      low = [spot.value-1000, 1].max
      high = spot.value+1000
      calls.each(spot.chr, low, high) do |c|
        if spot.watson?
          if c.position > spot.value and (p1.nil? or c.position < p1.position)
            p1 = c
          elsif c.position < spot.value and (m1.nil? or c.position > m1.position)
            m1 = c
          end
        else
          if c.position < spot.value and (p1.nil? or c.position > p1.position)
            p1 = c
          elsif c.position > spot.value and (m1.nil? or c.position < m1.position)
            m1 = c
          end
        end
      end
      
      if p1.nil? or m1.nil?
        invalid_coordinates += 1
      else
        f.puts "#{spot.chr}\t#{m1}\t#{p1}\t#{spot.id}\t#{spot.value}\t#{spot.strand}"
      end
    end
  end
end

puts "Skipped #{skipped} for invalid chromosome(s)" if skipped > 0
puts "Skipped #{invalid_coordinates} coordinate(s)" if invalid_coordinates > 0
