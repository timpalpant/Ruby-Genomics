require 'bio'
require 'genomic_data'

##
# Load an entire genome from a set of text files in a directory
# Loads the sequences as Bio::Sequence::NA
##
class Genome < GenomicData
  # Default genome assemblies for model organisms
  YEAST_GENOME = 'sacCer2'
  WORM_GENOME = 'notAvailable'
	FLY_GENOME = 'dm3'
  HUMAN_GENOME = 'notAvailable'
	
	# Resources directory
	RESOURCES = File.expand_path(File.dirname(__FILE__) + '/../resources')
	
	# Genome fasta directory
	GENOME_DIR = RESOURCES + '/genomes'
	  
  attr_reader :assembly
	
	# List the available genomes
	def self.list
		Dir.glob(GENOME_DIR + '/*').map { |g| File.basename(g)[0..-8] }
	end
	
	# Return the default yeast genome (sacCer2)
	def self.yeast
		self.load(YEAST_GENOME)
	end
	
  # Return the default worm genome
  def self.worm
    raise "The C. elegans (worm) genome is not available"
    #self.load(WORM_GENOME)
  end
	
	# Return the default fly genome
	def self.fly
		self.load(FLY_GENOME)
	end
  
  # Return the default human genome
  def self.human
    raise "The H. sapiens (human) genome is not available"
    #self.load(HUMAN_GENOME)
  end
	
  
	# Load a genome specified by name
	def self.load(name)
	  full_name = GENOME_DIR + '/' + name + '.genome'
	  raise "Specified genome assembly does not exist in resources/genomes/*" unless File.exist?(full_name)
	  
	  genome = nil
    File.open(full_name) do |f| 
      genome = Marshal.load(f)
    end
    
    return genome
	end
    
	
  # Create a new genome from a sequentially-numbered directory of Fasta files
  def initialize(dir)
		Dir.glob(GENOME_DIR + '/' + dir + '/*').each do |filename|
	    puts "Initializing #{dir} genome"
			
			if File.file?(filename)
    		# Load the genome data from fasta text files
	      lines = File.readlines(filename)
	      #header = lines.first.chomp
	      #chr = header[4..-1]
	      
				# Get the chromosome from the filename
				chr = File.basename(filename, '.raw')
				
				# Store as a BioRuby Nucleic Acid Sequence
	      self[chr] = Bio::Sequence::NA.new(lines.join)
        puts "#{chr}: #{filename} (#{self[chr].length} bases loaded)"
	    end
		end
		
		puts "Loaded #{self.length} chromosomes for assembly #{dir}"
	end
	
	# Return the number of base pairs in the genome
	# length() will return the number of chromosomes
	def bases
		self.values.inject(0) { |length,chr| length + chr.length }
	end

  def to_s
    str = "Genome #{@assembly.upcase}: containing #{self.bases} base pairs\n"
    self.each do |chr|
      str += "\tChromosome #{chr} (length: #{chr.length})\n"
    end
    
    return str
  end
  
  def composition
  	{'a'=>self.a_count, 't'=>self.t_count, 'g'=>self.g_count, 'c'=>self.c_count}
  end
	
  def a_count
  	count('a')
  end
  
	def a_content
    base_content('a')
	end
	
	def a_percent
	  (100 * a_content).to_f
	end
	
	def t_count
		count('t')
	end
	
  def t_content
   base_content('t')
  end
  
  def t_percent
    (100 * t_content).to_f
  end
 
  def g_count
  	count('g')
  end
  
  def g_content
   base_content('g')
  end
  
  def g_percent
    (100 * g_content).to_f
  end
 
  def c_count
  	count('c')
  end
  	
  def c_content
   base_content('c')
  end
  
  def c_percent
    (100 * c_content).to_f
  end
  
  
  private
  
  def base_content(base)
    Rational(count(base), self.bases)
  end
  
  def count(base)
		self.values.inject(0) { |count,chr| count + chr.count(base) }
  end
end