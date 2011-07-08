#!/usr/bin/env ruby1.9

# == Synopsis 
#   Find the average value of sequencing data around TTS, within genes, etc.
#
# == Usage 
#   Find the average value for list of chromosomal windows
#
#   regionIntensity.rb -i readcount.bw -w windows.bed -o RegionIntensities.txt
#
#   For help use: regionIntensity.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file(s) to average values (BigWig)
#   -w, --windows       List of windows to average around (in Bed format: chrXII  10345  10600)
#   -s, --statistic     Statistic to compute
#   -o, --output        Output file (flat list)
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
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i readcount.wig -w windows.bed -o WindowAverage.txt"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-w', '--window FILE', :required, "List of windows to average (in Bed format)" ) { |f| options[:windows] = f }
  opts.on( '-o', '--output FILE', :required, "Output file" ) { |f| options[:output] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches? or ARGV.length < 1
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Load the BigWig files
puts "\nInitializing BigWig file(s)" if ENV['DEBUG']
wigs = ARGV.map { |inputfile| WigFile.autodetect(inputfile) }

puts "\nComputing mean intensity for each region" if ENV['DEBUG']
File.open(options[:output], 'w') do |f|
  basenames = ARGV.map { |inputfile| File.basename(inputfile) }
  f.puts "#File\tStart-500..Start\tStart..Start+500\tStart-250..Start+250\tStart..Stop\tStart..Stop+500\tStart-500..Stop\tStart-500..Stop+500\tStop-500..Stop\tStop..Stop+500\tStop-250..Stop+250"
  
  totals = Array.new(wigs.length) { Array.new(10, 0) }
  count = 0
  BedFile.foreach(options[:windows]) do |spot|
    unless wigs.map { |wig| wig.include?(spot.chr) }.all?
      puts "Skipping spot because Wig(s) does not have data for #{chr}" if ENV['DEBUG']
      next
    end
    
    low = spot.watson? ? spot.start-500 : spot.start+500
    high = spot.watson? ? spot.stop+500 : spot.stop-500
    wigs.each_with_index do |wig,i|
      data = wig.query(spot.chr, low, high)
      
      # Start-500 .. Start
      totals[i][0] += data[0..500].mean
      
      # Start .. Start+500
      totals[i][1] += data[500..1000].mean
      
      # Start-250 .. Start+250
      totals[i][2] += data[250..750].mean
      
      # Start .. Stop
      totals[i][3] += data[500..-500].mean
      
      # Start .. Stop+500
      totals[i][4] += data[500..-1].mean
      
      # Start-500 .. Stop
      totals[i][5] += data[0..-500].mean
      
      # Start-500 .. Stop+500
      totals[i][6] += data.mean
      
      # Stop-500 .. Stop
      totals[i][7] += data[0..-500].mean
      
      # Stop .. Stop+500
      totals[i][8] += data[-500..-1].mean
      
      # Stop-250 .. Stop+250
      totals[i][9] += data[-750..-250].mean
      
    end
    
    count += 1
  end
  
  avgs = totals.map { |w| w.map { |t| t / count } }
  
  basenames.each_with_index do |name,i|
    f.puts "#{name}\t" + avgs[i].join("\t")
  end
end
