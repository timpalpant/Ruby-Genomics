#!/usr/bin/env ruby1.9

# == Synopsis 
#   Compare called nucleosomes around specific loci
#
# == Usage 
#  Compare nucleosomes from -3 to 3 around the loci listed in orfs.bed
#
#   overlapFinder.rb -1 wt-nukes.txt -2 mutant-nukes.txt -o overlappingNukes.txt
#
#   For help use: overlapFinder.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -1, --wildtype      Nuke calls from wildtype / sample 1
#   -2, --mutant        Nuke calls from mutant / sample 2
#   -b, --bases         Number of bases required for overlap
#   -o, --output        Nuke comparison results
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
require 'single_bp_data'
require 'nucleosome'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -1 wt-nukes.txt -2 mutant-nukes.txt -o nuke-comparison.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-1', '--wildtype FILE', :required, "Nuke calls from sample 1" ) { |f| options[:wt] = f }
  opts.on( '-2', '--mutant FILE', :required, "Nuke calls from sample 2" ) { |f| options[:mutant] = f }
  options[:overlap] = 73
  opts.on( '-b', '--bases N', "Number of bases required for overlap (default: 73)" ) { |n| options[:overlap] = n.to_i }
  opts.on( '-o', '--output FILE', :required, "Output summary file" ) { |f| options[:output] = f }
      
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

puts 'Loading nucleosomes from wildtype' if ENV['DEBUG']
wildtype = NukeCalls.load(options[:wt])
puts 'Loading nucleosomes from mutant' if ENV['DEBUG']
mutant = NukeCalls.load(options[:mutant])
  
overlap = GenomicData.new
mutant_chr, current_chr = nil, nil
wildtype.each do |chr,nukes|
	puts "Processing nucleosomes on chromosome #{chr}" if ENV['DEBUG']
	overlap[chr] ||= Array.new
	
	nukes.each do |wt_nuke|
    # Look for an overlapping nuke in the mutant
    overlapping_nukes = mutant_chr.select do |mut_nuke|
      if mut_nuke.start > wt_nuke.start
        (wt_nuke.stop - mut_nuke.start) > options[:overlap]
      else
        (mut_nuke.stop - wt_nuke.start) > options[:overlap]
      end
    end
    
    best_match = nil
    if overlapping_nukes.length > 1
      puts "Warning! More than one nucleosome overlaps by #{options[:overlap]} base pairs!"
      best_match = overlapping_nukes.sort { |nuke1, nuke2|
        nuke1_overlap = if nuke1.start > wt_nuke.start
          wt_nuke.stop - nuke1.start
        else
          nuke1.stop - wt_nuke.start
        end
        
        nuke2_overlap = if nuke2.start > wt_nuke.start
          wt_nuke.stop - nuke2.start
        else
          nuke2.stop - wt_nuke.start
        end
        
        nuke1_overlap <=> nuke2_overlap
      }.reverse.first
    elsif overlapping_nukes.length == 1
      best_match = overlapping_nukes.first
    end
    
    unless best_match.nil?
      overlap[chr] << [wt_nuke, best_match]
    end
	end
end
puts "Found #{overlap.length} overlapping nucleosomes" if ENV['DEBUG']

puts 'Writing overlapping nucleosomes to output file'
File.open(options[:output],'w') do |out|
  out.puts "chrom.A\tStart.A\tStop.A\tDyad.A\tDyad.Stdev.A\tDyad.Count.A\tchrom.B\tStart.B\tStop.B\tDyad.B\tDyad.Stdev.B\tDyad.Count.B\tOverlap"
  overlap.each do |nuke_pair|
    wt_nuke = nuke_pair.first
    mut_nuke = nuke_pair.last
    
    overlap = if mut_nuke.start > wt_nuke.start
      wt_nuke.stop - mut_nuke.start
    else
      mut_nuke.stop - wt_nuke.start
    end
    
    out.puts "{wt_nuke.chr}\t#{wt_nuke.start}\t#{wt_nuke.stop}\t#{wt_nuke.dyad}\t#{wt_nuke.dyad_stdev}\t#{wt_nuke.dyad_count}\t{mut_nuke.chr}\t#{mut_nuke.start}\t#{mut_nuke.stop}\t#{mut_nuke.dyad}\t#{mut_nuke.dyad_stdev}\t#{mut_nuke.dyad_count}\t#{overlap}"
  end
end
  
global_stop = Time.now
puts "Time elapsed: #{global_stop-global_start}" if ENV['DEBUG']