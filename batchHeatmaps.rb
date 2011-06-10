#!/usr/bin/env ruby1.9

# == Synopsis 
#   Generate heatmaps for the product of a set of sequencing data and alignment coordinates
#
# == Usage 
#   Generate heatmaps for the data in seqdir against the alignment points in aligns
#
#   batchHeatmaps.rb -d seqs -l aligns -o output_dir
#
#   For help use: wigCorrelate.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -d, --data          Sequencing data
#   -l, --loci          Loci to align to
#   -o, --output        Output directory
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
require 'forkmanager'
require 'unix_file_utils'


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
  opts.on( '-d', '--data DIR', :required, "Sequencing data directory" ) { |f| options[:data] = File.expand_path(f) }
  opts.on( '-l', '--loci DIR', :required, "Directory of alignment files" ) { |f| options[:loci] = File.expand_path(f) }
  opts.on( '-k', '--keep-matrices', "Keep aligned heatmap matrices" ) { |b| options[:keep] = b }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of threads to use (default: 2)" ) { |n| options[:threads] = n.to_i }
  options[:range] = '-3:3'
  opts.on( '-r', '--range MIN:MAX', "Range to use for heatmaps (default: -3:3)" ) { |s| options[:range] = s }
  opts.on( '-o', '--output DIR', "Output directory (default: loci directory)" ) { |f| options[:output] = File.expand_path(f) }
  
	# Parse the command-line arguments
	opts.parse!
	
	# Validate the required parameters
	if opts.missing_switches?
	  puts opts.missing_switches
	  puts opts
	  exit
	end
  
  # Default output directory is the directory with the alignment loci
  options[:output] = options[:loci] if options[:output].nil?
  options[:output] = File.dirname(options[:output]) if File.file?(options[:output])  

  # Check that matrix2png is in the PATH
  raise "Cannot find executable: matrix2png in $PATH which is required to make heatmaps" if File.which('matrix2png').nil?
end


##
# SETUP
##

# Parse the range to use for heatmap colors
range = options[:range].split(':')
raise "Invalid range! Should be of the form MIN:MAX" if range.length != 2
min = range.first.to_f
max = range.last.to_f
raise "Invalid range! Should be of the form MIN:MAX" if min >= max

# Queue the sequencing data files to process
data_files = if File.file?(options[:data])
  [options[:data]]
elsif File.directory?(options[:data])
  Dir.glob(options[:data]+'/*').select { |f| File.file?(f) and File.extname(f) == '.bw' or File.extname(f) == '.bigwig' }
else
  raise "Sequencing data parameter is neither a single file nor a directory! Exiting..."
end
puts "\nQueued #{data_files.length} sequencing data file(s) to process:"
data_files.each_with_index { |f,i| puts "#{i+1}\t#{File.basename(f)}" }

# Queue the alignment files to process
loci_files = if File.file?(options[:loci])
  [options[:loci]]
elsif File.directory?(options[:loci])
  Dir.glob(options[:loci]+'/*').select { |f| File.file?(f) and File.extname(f) == '.bed' }
else
  raise "Alignment loci parameter is neither a single file nor a directory! Exiting..."
end
puts "\nQueued #{loci_files.length} alignment loci file(s) to process:"
loci_files.each_with_index { |f,i| puts "#{i+1}\t#{File.basename(f)}" }

# Compute how many heatmaps we will be generating
puts "\nWill be generating #{data_files.length * loci_files.length} heatmaps"

# Initialize the forkmanager for running in parallel
pm = Parallel::ForkManager.new(options[:threads])


##
# HEATMAP GENERATION
##

loci_files.each do |loci|
  puts "\nProcessing alignment loci #{File.basename(loci)}"
  
  # Make a directory for the output files with the same name as the loci file
  heatmap_dir = options[:output] + '/' + File.basename(loci, '.bed')
  Dir.mkdir(heatmap_dir) unless File.directory?(heatmap_dir)
  
  # Generate a heatmap for each sequencing dataset
  data_files.each do |input|
    # Run each heatmap job in a separate process
    pm.start and next
    
    puts "\tdata file: #{File.basename(input)}"
    
    # Align values in a matrix for a heatmap
    aligned_matrix = heatmap_dir + '/' + File.basename(input) + '.matrix'
    %x[ ruby1.9 visualization/matrixAligner.rb -i #{input} -l #{loci} -o #{aligned_matrix} ]
    
    # Generate a heatmap with matrix2png (must be on the PATH)
    heatmap = heatmap_dir + '/' + File.basename(File.basename(input, '.bw'), '.bigwig') + '.png'
    %x[ matrix2png -data #{aligned_matrix} -range #{min}:#{max} -mincolor blue -maxcolor yellow -bkgcolor white -missingcolor gray -title "#{File.basename(input)} aligned to #{File.basename(loci)} (+/- 500 bp)" -b -s > #{heatmap} ]
    
    # Remove the aligned matrix unless we are keeping them
    File.delete(aligned_matrix) unless options[:keep]

    pm.finish(0)
  end
end

# Wait for all of the heatmap jobs to finish
pm.wait_all_children

puts "Complete!"
