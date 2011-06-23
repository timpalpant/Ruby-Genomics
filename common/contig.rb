require 'fixed_precision'
require 'genomic_index_error'
require 'unix_file_utils'

##
# Represents a contiguous block of genomic data values
# Returned when querying a WigFile / SpotArray
#
# Since chromosomal coordinates can always be indicated by integers, store as an Array
##
class Contig < Array
  attr_accessor :chr, :start, :step, :span

  def initialize(length = 0, chr = 'unknown', start = 1, step = 1, span = 1)
    super(length, 0)
    @chr = chr
    @start = start
    @step = step
    @span = span
  end
  
  # Parse the header arguments of a fixedStep/variableStep line
  def self.parse_wig_header(line)
    unless line.chomp.start_with?('fixedStep') or line.chomp.start_with?('variableStep')
      raise ContigError, "Not a valid Wig Contig header"
    end
  
    chrom, start, step, span = 'unknown', 1, 1, 1
    line.chomp.split(' ').each do |opt|
      keypair = opt.split('=')
      key = keypair.first
      value = keypair.last
      
      case key
        when 'chrom'
          chrom = value
        when 'start'
          start = value.to_i
        when 'step'
          step = value.to_i
        when 'span'
          span = value.to_i
      end
    end

    return self.new(0, chrom, start, step, span)
  end
  
  # Load a chromosome from a specific section of a Wig file
  # (first line should be a fixedStep/variableStep line)
  def self.load_wig(data_file, start_line, end_line)
    # Grab all the lines
    data = File.lines(data_file, start_line, end_line).reject { |line| line.chomp.empty? }
      
    # Parse the fixedStep/variableStep header
    header = data.first
    chr = parse_wig_header(header)

    # Parse all the values into floats
    if header.start_with?('variableStep')
      # Reconfigure the data if it was variableStep
      chr.start = data[1].split("\t").first.to_i
      chr.step = 1
      chr.span = 1
    
      prev_base, prev_value = nil, nil
      data[1..-1].each do |line|
        entry = line.chomp.split("\t")
        base = entry.first.to_i
        value = entry.last.to_f
        
        if prev_base and prev_value
          for bp in prev_base-chr.start..base-chr.start
            chr[bp] = prev_value
          end
        end
        
        prev_base = base
        prev_value = value
      end
    elsif header.start_with?('fixedStep')
      chr[0..-1] = data[1..-1].map { |line| line.to_f }
    else
      raise "Wig chromosome header is neither fixedStep nor variableStep!"
    end
    
    return chr  
  end

  # The last base pair with data
  def stop
    @start + self.length - 1
  end
  
  # If this Chromosome contains data for the specified range
  def include?(low, high = self.length-1)
    low >= start and low <= stop and high >= start and high <= stop
  end
  
  # Get a subsequence of data
  def bases(from, to)
    low_bp = [from, to].min
    high_bp = [from, to].max
    
    raise ContigError, "Chromosome does not include bases #{low_bp}..#{high_bp}" unless self.include?(low_bp, high_bp)
    
    data = self[low_bp-@start..high_bp-@start]
    data = data.reverse if from > to
    return data
  end
  
  # Alias for bases
  def subseq(low, high)
    bases(low, high)
  end
  
  # Output this Contig as a fixedStep list of values (one per line)
  def to_s
    str = StringIO.new("fixedStep chrom=#{@chr}")
    
    str << " start=#{@start}" if @start
    str << " step=#{@step}" if @step
    str << " span=#{@span}" if @span

    self.each { |value| str << "\n" << value.to_s(5) }
    
    return str.string
  end
end


# For converting an Array to a Chromosome
class Array
  def to_contig(chr = 'unknown', start = 1, step = 1, span = 1)
    chr = Contig.new(self.length, chr, start, step, span)
    chr[0..-1] = self
    return chr
  end
end


# Raised if something goes wrong with a Chromosome
class ContigError < StandardError
end
