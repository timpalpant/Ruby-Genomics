#
#  Assembly.rb
#  BioRuby
#
#  Created by Timothy Palpant on 5/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'genomic_data'

class Assembly < GenomicData
  # Default assemblies for model organisms
  YEAST_GENOME = 'sacCer2'
  WORM_GENOME = 'ce9'
  FLY_GENOME = 'dm3'
  HUMAN_GENOME = 'hg19'
  
  # Resources directory
  RESOURCES = File.expand_path(File.dirname(__FILE__) + '/../resources')
  
  # Genome fasta directory
  ASSEMBLY_DIR = RESOURCES + '/assemblies'
  
  attr_reader :name, :len_file
  
  # Initialize a new assembly
  def initialize(name, len_file = nil)
    @name = name
    @len_file = len_file
  end
  
  # Write this assembly file to disk in len format (chr_id \t num_bases)
  def to_len(filename)
    File.open(File.expand_path(filename), 'w') do |f|
      self.each do |chr, chr_len|
        f.puts "#{chr}\t#{chr_len}"
      end
    end
  end
  
  # Summary info about this Assembly
  def to_s
    str = "Assembly: #{@name}"
    self.each do |chr, chr_length|
      str += "\n\t#{chr}: #{chr_length} bases"
    end
    
    return str
  end

  # List the available assemblies
  def self.list
    Dir.glob(ASSEMBLY_DIR + '/*').map { |g| File.basename(g)[0..-5] }
  end

  # Return the default yeast genome (sacCer2)
  def self.yeast
    self.load(YEAST_GENOME)
  end

  # Return the default worm genome
  def self.worm
    self.load(WORM_GENOME)
  end

  # Return the default fly genome
  def self.fly
    self.load(FLY_GENOME)
  end

  # Return the default human genome
  def self.human
    self.load(HUMAN_GENOME)
  end

  # Load an assembly specified by name
  def self.load(name)
    len_file = ASSEMBLY_DIR + '/' + name + '.len'
    raise "Specified genome assembly does not exist in resources/assemblies/*" unless File.exist?(len_file)

    assembly = self.new(name, len_file)
    File.foreach(len_file) do |line|
      entry = line.chomp.split("\t")
      raise "Invalid entry in #{full_name}" if entry.length != 2
      
      assembly[entry[0]] = entry[1].to_i
    end

    return assembly
  end
  
end