#!/usr/bin/env ruby1.9

# == Synopsis 
#   Correct for background signal in competition ChIP-chip experiments
#
# == Author
#   Timothy Palpant
#
# == Copyright
#   Copyright (c) 2011 Timothy Palpant. Licensed under the GPL v3:
#   http://www.gnu.org/licenses/gpl-3.0.txt

COMMON_DIR = File.expand_path(File.dirname(__FILE__) + '/common')
$LOAD_PATH << COMMON_DIR unless $LOAD_PATH.include?(COMMON_DIR)
require 'bundler/setup'
require 'pickled_optparse'
require 'bio-genomic-file'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -d seqs -l aligns -o output_dir"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-d', '--data DIR', :required, "Data directory" ) { |f| options[:data] = File.expand_path(f) }
  opts.on( '-b', '--background FILE1,FILE2,...', Array, :required, "Files to use for background estimation" ) { |a| options[:bgd] = a.map { |f| File.expand_path(f) } }
  options[:suffix] = String.new
  opts.on( '-s', '--suffix S', "Suffix for normalized files" ) { |s| options[:suffix] = s }
  opts.on( '-o', '--output DIR', "Output directory (default: data directory)" ) { |f| options[:output] = File.expand_path(f) }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
  
  if not File.directory?(options[:data])
    puts "Data directory is invalid!"
    exit(-1)
  end
  
  # By default, dump data in the same directory and include a suffix
  if not options.include?(:output)
    options[:output] = options[:data]
    options[:suffix] = '.bgdCorrected' if not options.include?(:suffix)
  end
  
  if not File.directory?(options[:output])
    puts "Output directory is invalid!"
    exit(-1)
  end
end

# Estimate the background
total = Hash.new(0)
count = Hash.new(0)
options[:bgd].each do |f|
  BedGraphFile.foreach(f) do |entry|
    bgd[entry.query_string] += entry.value 
    count[entry.query_string] += 1
  end
end

# Compute the average for each spot
avg = Hash.new
total.each_key { |k| avg[k] = total[k] / count[k] }

# Compute the background constant for each spot
c = Hash.new
avg.each_key do |k| 
  ratio = 2**avg[k]
  c[k] = ratio / (1 - ratio)
end
puts "Computed the background from #{options[:bgd]} time points for #{c.length} spots" if ENV['DEBUG']

# Normalize each of the time points to the background
Dir.glob(options[:data] + '/*').each do |f|
  ext = File.extname(f)
  basename = File.basename(f, ext)
  File.open(options[:output] + "/#{basename}#{options[:suffix]}#{ext}", 'w') do |output|
    # Write the BedGraph header
    output.puts Utils::UCSC::TrackHeader.new(:type => 'bedGraph', 
                                             :name => basename, 
                                             :description => basename)
  
    BedGraphFile.foreach(f) do |entry|
      k = entry.query_string
      m = 2**entry.value
    
      # Correct the measured value for background and re-log2 transform
      corrected_ratio = (c[k]*(1 - m) - m) / (c[k]*(m - 1) - 1)
      entry.value = Math.log(corrected_ratio, 2)
      
      # Write the corrected spot to disk
      output.puts entry.to_bedgraph
    end
  end
end
