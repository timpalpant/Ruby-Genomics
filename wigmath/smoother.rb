#!/usr/bin/env ruby1.9

# == Synopsis 
#   Gaussian smooths (Big)Wig Files using FFT-based convolution
#
# == Usage 
#   Smooth file1.bw:
#
#   smoother.rb -i file1.bw -o file1.smoothed.bw
#
#   For help use: smoother.rb -h
#
# == Options
#   -h, --help          Displays help message
#   -i, --input         Input file
#   -s, --stdev         Standard deviation of the Gaussian
#   -w, --window        Window size
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
require 'reference_assembly'
require 'wig_transform'
require 'fftw3'
include Bio

# This hash will hold all of the options parsed from the command-line by OptionParser.
options = Hash.new
ARGV.options do |opts|
  opts.banner = "Usage: ruby #{__FILE__} -i file1.bw -o file1.smoothed.wig"
  # This displays the help screen, all programs are assumed to have this option.
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts opts
    exit
  end
  
  # List all parameters
  opts.on( '-i', '--input FILE', :required, "Input (Big)Wig file" ) { |f| options[:input] = f }
  options[:sdev] = 20
  opts.on( '-s', '--sdev NUM', "Standard deviation of the Gaussian in base pairs (default 20)" ) { |num| options[:sdev] = num.to_i }
  options[:window_size] = 3
  opts.on( '-w', '--window NUM', "Number of standard deviations +/- to make a window (default 3)" ) { |num| options[:window_size] = num.to_i }
  options[:threads] = 2
  opts.on( '-p', '--threads N', "Number of processes (default: 2)" ) { |n| options[:threads] = n.to_i }
  opts.on( '-g', '--genome ASSEMBLY', :required, "Genome assembly" ) { |g| options[:genome] = g }
  opts.on( '-o', '--output FILE', :required, "Output file (BigWig)" ) { |f| options[:output] = f }
      
  # Parse the command-line arguments
  opts.parse!
  
  # Validate the required parameters
  if opts.missing_switches?
    puts opts.missing_switches
    puts opts
    exit
  end
end

# Smoothing requires padding on either end
padding = options[:sdev] * options[:window_size]

# Initialize the wig files to smooth
wig = WigFile.autodetect(options[:input])

# Initialize the output assembly
assembly = ReferenceAssembly.load(options[:genome])

##
# Run the smoothing transformation
##
filter = nil
wig.transform(options[:output], assembly, :in_processes => options[:threads]) do |chr, chunk_start, chunk_stop|
  # Don't pad off the end of the chromosome
  query_start = [wig.chr_start(chr), chunk_start-padding].max
  query_stop = [chunk_stop+padding, wig.chr_stop(chr)].min
  
  # Actual padding
  padding_left = chunk_start - query_start
  padding_right = query_stop - chunk_stop
  
  # Query for the values from the Wig file
  chunk = wig.query(chr, query_start, query_stop)
  
  # Construct the filter if necessary (cache for performance)
  # TODO: Construct the filter in frequency space directly
  # since the Fourier transform of a Gaussian is another Gaussian
  if filter.nil? or filter.length != chunk.length
    puts "Generating Gaussian filter with length #{chunk.length}" if ENV['DEBUG']
    raise "Cannot smooth a chunk with size #{chunk.length} with a filter of size #{2*padding+1}!" if 2*padding+1 > chunk.length
    
    gaussian = NArray.float(chunk.length)
    for x in -padding..padding
      gaussian[gaussian.length/2 + x] = Math.exp(-((x**2)/(2*(options[:sdev]**2))))
    end
    gaussian /= gaussian.sum
    filter = FFTW3.fft(gaussian)
  end
  
  # Convolve the filter with the data by multiplying in frequency space and then inverting the Fourier transform
  p = NArray.float(chunk.length)
  chunk.each do |bp, value|
    p[bp-chunk.start] = value
  end
  f = FFTW3.fft(p)
  s = FFTW3.ifft(filter * f).real / f.length
  
  # ifftshift
  h = s.length / 2
  for i in 0...h
    tmp = s[i]
    s[i] = s[h+i]
    s[h+i] = tmp
  end
  
  # Store the convolved results in a new Contig
  smoothed = Genomics::Contig.new(chr)
  for i in padding_left...s.length-padding_right
    smoothed.set(chunk.start+i, s[i])
  end
  
  smoothed
end
