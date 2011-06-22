#!/usr/bin/env ruby1.9

# == Synopsis 
#   Marks peak summits relative to features
#
# == Usage 
#   Mark peak summits on a scatter plot relative to TSS
#
#   markBindingSites.rb -i Nhp6A-summits.bed -l TSS-500.bed -o Nhp6A-TSS.peaks.txt
#
#   For help use: markBindingSites.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file with peak summits
#   -l, --loci          List of windows to align (chr  start  stop  id   alignmentBase)
#   -o, --output        Output file with peak summits relative to loci alignment point
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
require 'macs'

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i Nhp6A-summits.bed -l TSS-500.bed -o Nhp6A.peaks.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input file with peak calls (MACS output)" ) { |f| options[:input] = f }
  opts.on( '-l', '--loci FILE', :required, "List of loci to align to (Bed format)" ) { |f| options[:loci] = f }
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


puts "Loading the list of alignment loci" if ENV['DEBUG']
loci = BedFile.load(options[:loci])
puts "#{loci.num_spots} alignment loci"
  
puts "\nLoading binding sites" if ENV['DEBUG']
bs = MACSFile.load(options[:input])
puts "#{bs.num_spots} binding sites"

puts "\nMarking binding sites relative to alignment loci" if ENV['DEBUG']
File.open(options[:output], 'w') do |f|
  f.puts "#Chromosome\tStart\tStop\tAlignment Point\tBinding Site\tRelative Binding Site\tTags"
  
  bs.each do |chr,sites|
    if not loci.include?(chr)
      puts "Skipping binding sites because no loci for #{chr}"
      next
    end
    
    sites.each do |site|
      selected_loci = loci[chr].select { |spot| spot.include?(site.summit) }
  
      if selected_loci.length > 0
        selected_loci.each do |locus|
          offset = if locus.watson?
            site.summit - locus.value.to_i
          else
            locus.value.to_i - site.summit
          end
          
          f.puts "#{chr}\t#{locus.start}\t#{locus.stop}\t#{locus.value.to_i}\t#{site.summit}\t#{offset}\t#{site.tags}"
        end
      else
        puts "No spots found for binding site (#{chr}#{site})" if ENV['DEBUG']
      end
    end
  end
end
