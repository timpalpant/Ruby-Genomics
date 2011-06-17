require 'bio'
require 'genomic_data'
require 'assembly'

##
# Load an entire genome from a set of Fasta files in a directory
# Loads the sequences as Bio::Sequence::NA
##
class Genome < GenomicData	
	# Load a genomic reference sequence from a fasta file, or a directory of fasta files
  def initialize(filename)
    expanded = File.expand_path(filename)
    if File.file?(expanded)
      Bio::FlatFile.foreach(Bio::FastaFormat, expanded) do |entry|
        self[entry.definition] = entry.naseq
      end
    elsif File.directory?(expanded)
      Dir.glob(expanded + '/*').each do |subfile|        
        if File.file?(subfile)
          Bio::FlatFile.foreach(Bio::FastaFormat, subfile) do |entry|
            self[entry.definition] = entry.naseq
          end
        end
      end
    else
      raise "Genome to load must be either a single Fasta file or a directory of Fasta files, one per chromosome"
    end
		
		puts "Loaded #{self.length} chromosomes (#{self.bases} base pairs) from #{File.basename(filename)}" if ENV['DEBUG']
	end
	
	# Return the number of base pairs in the genome
	# length() will return the number of chromosomes
	def bases
		self.values.inject(0) { |length,chr| length + chr.length }
	end
  
  # Reduce this genome to an Assembly object (just chromosome id's and their lengths)
  def to_assembly(name)
    a = Assembly.new(name, nil)
    
    self.each do |chr, seq|
      a[chr] = seq.length
    end
    
    return a
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