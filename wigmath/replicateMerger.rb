#!/usr/bin/env ruby1.9

# == Synopsis 
#   Normalize a Wig file to input and average multiple replicates
#
# == Examples
#   This command processes seqData.bw:
#     replicateMerger.rb -g sacCer2 -r output.relative.wig -a output.absolute.wig \
#       file1.IP.wig,file1.input.wig file2.IP.wig,file2.input.wig [...]
#
#   For help use: replicateMerger.rb -h
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
require 'bio-genomic-file'
require 'pickled_optparse'
require 'reference_assembly'
require 'wig_transform'
require 'utils/parallelizer'
require 'fileutils'
require 'rbconfig'
include Bio

RUBY_INTERPRETER = File.join(Config::CONFIG["bindir"],
                             Config::CONFIG["RUBY_INSTALL_NAME"] +
                             Config::CONFIG["EXEEXT"])
DIFFERENCE_SCRIPT = File.expand_path(File.dirname(__FILE__) + '/difference.rb')
DIVIDE_SCRIPT = File.expand_path(File.dirname(__FILE__) + '/divide.rb')
LOGGER_SCRIPT = File.expand_path(File.dirname(__FILE__) + '/logger.rb')
ZSCORER_SCRIPT = File.expand_path(File.dirname(__FILE__) + '/zscorer.rb')
AVERAGER_SCRIPT = File.expand_path(File.dirname(__FILE__) + '/averager.rb')

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -a absolute.wig -r relative.wig input1.IP.wig,input1.input.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # Input/output arguments
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-g', '--genome ASSEMBLY', :required, "Genome assembly" ) { |g| options[:genome] = g }
  options[:base] = 2
  opts.on( '-b', '--base N', "Logarithm base (default: 2)" ) { |n| options[:base] = n.to_i }
  opts.on( '-a', '--absolute FILE', :required, "Absolute difference output Wig file" ) { |f| options[:absolute] = f }
  opts.on( '-r', '--relative FILE', :required, "Relative difference output Wig file" ) { |f| options[:relative] = f }
  
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Process the input files
DataPair = Struct.new(:signal, :input)
replicates = ARGV.map do |arg| 
  entry = arg.split(',')
  raise "Input files not specified correctly!" if entry.length != 2
  DataPair.new(entry[0], entry[1])
end

# Ensure that all intermediate files are deleted if something crashes
NormalizedPair = Struct.new(:absolute, :relative)
tmp_files = Array.new
begin
  i = 0
  processed = replicates.p_map(:in_threads => options[:threads])  do |pair|
    i += 1    

    # Compute the absolute difference and Z-score
    tmp_difference = "replicate#{i}.diff.wig"
    tmp_files << tmp_difference
    %x[ #{RUBY_INTERPRETER} #{DIFFERENCE_SCRIPT} -p #{options[:threads]} -g #{options[:genome]} -m #{pair.signal} -s #{pair.input} -o #{tmp_difference}]
    absolute = tmp_difference
    #absolute = "replicate#{i}.diff.zscored.wig"
    #tmp_files << absolute
    #%x[ #{RUBY_INTERPRETER} #{ZSCORER_SCRIPT} -p #{options[:threads]} -g #{options[:genome]} -i #{tmp_difference} -o #{absolute} ]
  
    # Compute the relative difference, log-transform and Z-score
    tmp_divide = "replicate#{i}.div.wig"
    tmp_files << tmp_divide
    %x[ #{RUBY_INTERPRETER} #{DIVIDE_SCRIPT} -p #{options[:threads]} -g #{options[:genome]} -1 #{pair.signal} -2 #{pair.input} -o #{tmp_divide}]
    tmp_log = "replicate#{i}.div.log.wig"
    tmp_files << tmp_log
    %x[ #{RUBY_INTERPRETER} #{LOGGER_SCRIPT} -p #{options[:threads]} -g #{options[:genome]} -b #{options[:base]} -i #{tmp_divide} -o #{tmp_log}]
    relative = tmp_log
    #relative = "replicate#{i}.div.log.zscored.wig"
    #tmp_files << relative
    #{}%x[ #{RUBY_INTERPRETER} #{ZSCORER_SCRIPT} -p #{options[:threads]} -g #{options[:genome]} -i #{tmp_log} -o #{relative} ]
    
    NormalizedPair.new(absolute, relative)
  end

  # If more than one replicate provided, average replicates
  # Otherwise, just move the single replicate into the output files
  if replicates.length > 1
    # Average absolute differences
    absolute_replicates = processed.map { |p| p.absolute }.join(' ')
    %x[ #{RUBY_INTERPRETER} #{AVERAGER_SCRIPT} -p #{options[:threads]} -g #{options[:genome]} -o #{options[:absolute]} #{absolute_replicates} ]
    
    # Average relative differences
    relative_replicates = processed.map { |p| p.relative }.join(' ')
    %x[ #{RUBY_INTERPRETER} #{AVERAGER_SCRIPT} -p #{options[:threads]} -g #{options[:genome]} -o #{options[:relative]} #{relative_replicates} ]
  else
    FileUtils.move(processed[0].absolute, options[:absolute])
    FileUtils.move(processed[0].relative, options[:relative])
  end
ensure
  tmp_files.each do |f|
    File.delete(f) if File.exist?(f)
  end
end
