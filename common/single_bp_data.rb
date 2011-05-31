require 'genomic_data'
require 'chromosome'
require 'single_bp_math'

##
# Single-base pair resolution genomic data stored as Arrays
##
class SingleBPData < GenomicData
	include SingleBPMath
	
  attr_accessor :name, :description, :assembly
  
  # Return a new empty SingleBPData (zeros) with dimensions of an assembly
  def self.for_assembly(a)    
    empty_seq = self.new
		empty_seq.assembly = a.name
  	a.each { |chr,num_bases| empty_seq[chr] = Chromosome.new(num_bases) }
    return empty_seq
  end
  
  # Summary information about this sequencing file object
  def to_s
    out = "Sequencing File: #{@name}\n"
    out << "Assembly: #{@assembly}\n"
    out << "Total number of values: #{self.num_values}\n"
    out << "Total: #{self.total}\n"
    out << "Mean: #{self.mean}\n"
    out << "StDev: #{self.stdev}\n"
    
    # Summary for each Chromosome
  	out << self.map { |chr_id,chr_data| "Chromosome #{chr_id}: #{chr_data.length} values" }.join("\n")
      
    return out
  end
  
  # Output the data in Wig format
  def to_wig(filename)
    File.open(filename,'w') do |f|
      # TODO: should be rewritten to intelligently use step size
      f.puts Wig.track_header(@name, @description) 
          
      self.each do |chr_id,values|
        f.puts Wig.fixed_step(chr_id, values)
        f.puts values
      end
    end
  end
end


##
# Load a fixedStep Wig file wholly into memory as a SingleBPData
##
class Wig < SingleBPData
  def self.load(filename)
    wig = self.new
    
		# Piggyback on WigFile parsing
    # Iterate over each chromosome and store permanently
    cached_wig = WigFile.new(filename)
    wig.name = cached_wig.name
    wig.description = cached_wig.description
    
    cached_wig.each do |chr,values|
    	wig[chr] = values
    end
    
    return wig
  end
  
  # Make a Wig track header with the given name and description
  def self.track_header(name = '', description = '')
  	"track type=wiggle_0 name=\"#{name}\" description=\"#{description}\" autoScale=\"off\" visibility=\"full\""
  end
  
  # Make a Wig fixedStep chromosome header for the given ID and Chromosome
  def self.fixed_step(chr_id, chr = nil)
  	str = "fixedStep chrom=#{chr_id}"
  	
  	# fixedStep parameters
  	unless chr.nil?
  	  str << " start=#{chr.start}" if chr.start
  	  str << " step=#{chr.step}" if chr.step
	    str << " span=#{chr.span}" if chr.span
  	end
  	
  	return str
  end
end

class WigError < StandardError
end