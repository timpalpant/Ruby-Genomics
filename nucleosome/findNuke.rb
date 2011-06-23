#!/usr/bin/env ruby1.9

# == Synopsis 
#   Finds the first nucleosome from either the 5' or 3' end of a window (such as an ORF)
#
# == Usage 
#   Finds the first nucleosome for the spots in orfs.bed for the nukes in nukes.txt
#
#   findNuke.rb -i nukes.txt -l orfs.bed -o orfs.nukes.bed
#
#   For help use: findNuke.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Nuke calls
#   -g, --genome        Genome assembly to use
#   -l, --loci          Loci to search for nukes in
#   -o, --output        End-nuke calls
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
require 'bed'
require 'nucleosome'

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
  options[:reverse] = false
  opts.on( '-r', '--reverse', "Search from the 3' ends of windows (default: false)" ) { |b| options[:reverse] = b }  
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


# Load the windows
skipped, invalid_coordinates = 0, 0
File.open(options[:output], 'w') do |f|
  BedFile.open(options[:loci]) do |bed|
    NukeCallsFile.open(options[:input]) do |calls|
      direction = options[:reverse] ? 3 : 5
      puts "Finding first nucleosome from #{direction}' end" if ENV['DEBUG']
      
      bed.each do |spot|
        if not calls.chromosomes.include?(spot.chr)
          skipped += 1
          next
        end
        
        first_nuke = nil
        calls.each(spot.chr, spot.start, spot.stop) do |nuke|
          # If this is the first nucleosome we've found, it automatically wins
          if first_nuke.nil?
            first_nuke = nuke.dyad
          # Otherwise it has to be the best
          else
            first_nuke = if (options[:reverse] and spot.watson?) or (not options[:reverse] and spot.crick?)
              nuke.dyad if nuke.dyad > first_nuke
            else
              nuke.dyad if nuke.dyad < first_nuke
            end
          end
        end
        
        if first_nuke == nil
          invalid_coordinates += 1
        else
          spot.value = first_nuke
          f.puts spot.to_bed
        end
      end
    end
  end
end

puts "No nucleosomes for #{skipped} chromosome(s)" if skipped > 0
puts "No nucleosomes for #{invalid_coordinates} coordinate(s)" if invalid_coordinates > 0
