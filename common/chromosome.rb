require 'gsl'
require 'fixed_precision'
require 'genomic_index_error'
require 'unix_file_utils'
require 'math_utils'
include GSL

##
# Represents a chromosome of single-bp genomic data
#   Returned when querying a WigFile for a specific chromosome with WigFile#chr(num)
#   or iterating over the chromosomes of a GenomicData
#
# Since chromosomal coordinates can always be indicated by integers, store as a GSL::Vector
##
class Chromosome < Vector
  attr_accessor :start, :step, :span

  def self.new(length, start = 1, step = 1, span = 1)
  	chr = self[length]
  	chr.start = start
  	chr.step = step
  	chr.span = span
  	return chr
  end
  
  # Load a chromosome from a specific section of a Wig file
  # (first line should be a fixedStep/variableStep line)
  def self.load(data_file, start_line, end_line)
		#num_lines = end_line - start_line + 1
		# If loading > 1 million values, split up and load in pieces
		# somewhat arbitrary, can be adjusted to accomodate memory footprint requirements
		#if num_lines > 1_000_000
			# Grab just the header line
			#header = File.lines(data_file, start_line, start_line).first
		#else
			# Grab all the lines
			data = File.lines(data_file, start_line, end_line).reject { |line| line.chomp.empty? }
			
			# Parse the fixedStep/variableStep header
			header = data.first
		
		# Allocate space for the data
		length = if header.start_with?('variableStep')
			data[-1].split("\t").first.to_i - data[1].split("\t").first.to_i + 1
		elsif header.start_with?('fixedStep')
			data[1..-1].length
		else
			raise ChromosomeError, "Chromosome header is not fixedStep or variableStep!"
		end
		chr = self[length]
		
		# Parse the header arguments
    header.split(' ').each do |opt|
      keypair = opt.split('=')
      key = keypair.first
      value = keypair.last
      
      case key
        when 'start'
          chr.start = value.to_i
        when 'step'
          chr.step = value.to_i
        when 'span'
          chr.span = value.to_i
      end
    end

    # Parse all the values into floats
		if header.start_with?('variableStep')
			# Reconfigure the data if it was variableStep
			chr.start = data[1].split("\t").first.to_i
			chr.step = 1
			chr.span = 1
		
			prev_base, prev_value = nil, nil
			data[1..-1].each do |line|
				entry = line.split("\t")
				base = entry.first.to_i
				value = entry.last.to_f
				
				if prev_base and prev_value
					chr[prev_base-chr.start..base-chr.start] = prev_value
				end
				
				prev_base = base
				prev_value = value
			end
		elsif header.start_with?('fixedStep')
			chr[0..-1] = data[1..-1].map { |line| line.to_f }
		end

		# Force GC
		data = nil
		GC.start
		
    return chr  
  end

	# The last base pair with data
  def stop
  	self.length + @start - 1
  end
  
  # If this Chromosome contains data for the specified range
  def include?(low, high = self.length-1)
		low >= start and low <= stop and high >= start and high <= stop
  end
  
  # Get a subsequence of data
  def bases(from, to)
		low_bp = Math.min(from, to)
		high_bp = Math.max(from, to)
		
		raise ChromosomeError, "Chromosome does not include bases #{low_bp}..#{high_bp}" unless self.include?(low_bp, high_bp)
		
  	data = self[low_bp-@start..high_bp-@start]
		data = data.reverse if from > to
		return data
  end
	
	# Alias for bases
	def subseq(low, high)
		bases(low, high)
	end
  
  # Output this chromosome as a fixedStep list of values (one per line)
  def to_s(particular = false)
  	# Weird artifacts sometimes crop up
  	return self.join("\n") unless particular
		
  	# If more control over the number of sig figs is needed (2-3x slower)
    self.to_a.map { |value| value.to_s(5) }.join("\n")
  end
end


# For converting an Array to a Chromosome
class Array
	def to_chr(start = 1, step = 1, span = 1)
		chr = Chromosome.alloc(self)
		chr.start = start
		chr.step = step
		chr.span = span
		return chr
	end
end


# Raised if something goes wrong with a Chromosome
class ChromosomeError < StandardError
end