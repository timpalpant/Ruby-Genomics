#
#  parallelizer.rb
#  BioRuby
#  Parallelize computations across chromosomes with ForkManager
#
#  Created by Timothy Palpant on 6/1/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'fileutils'
require 'tmpdir'
require 'forkmanager'
require 'ucsc_tools'

class Parallelizer
  def initialize(max_threads = 2)
    @pm = Parallel::ForkManager.new(max_threads, {'tempdir' => Dir.tmpdir})
  end
end

class WigComputationParallelizer < Parallelizer
  def initialize(output_file, chunk_size = 200_000, max_threads = 2)
    super(max_threads)
    @output = output_file
    @chunk_size = chunk_size
  end
  
  # Run a given block transform for each chromosome in chunks
  # The wig file provides the chromosomes and the chunk coordinates
  # Sort of a map-reduce approach
  # Return a WigFile handle to the output
  def run(wig)
    # Write the output file header
    header_file = @output+'.header'
    File.open(header_file, 'w') do |f|
      f.puts wig.track_header
    end

    # Keep track of all the temporary intermediate files (header first)
    tmp_files = [header_file]
  
    # Iterate over the Wig file chromosome-by-chromosome
    wig.chromosomes.each do |chr|
      chr_temp_file = @output+'.'+chr
      tmp_files << chr_temp_file      
      
      # Run in parallel processes managed by ForkManager
      @pm.start(chr) and next
      puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']

      # Write the chromosome fixedStep header
      File.open(chr_temp_file, 'w') do |f|
        f.puts Contig.new(0, chr, 1, 1, 1)
      end
      
      chunk_start = 1
      chr_length = wig.chr_length(chr)
      while chunk_start <= chr_length
        chunk_stop = [chunk_start+@chunk_size-1, chr_length].min
        puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
        
        output = yield(chr, chunk_start, chunk_stop)
        
        # Write this chunk to disk
        File.open(chr_temp_file, 'a') do |f|
          f.puts output.map { |value| value.to_s(5) }.join("\n")
        end
        
        chunk_start = chunk_stop + 1
      end

      @pm.finish(0)
    end
    
    # Wait for all of the child processes (each chromosome) to complete
    @pm.wait_all_children

    # Concatenate all of the temp file pieces into the final output
    File.cat(tmp_files, @output)

    # Delete the individual temp files created by each process
    tmp_files.each { |filename| File.delete(filename) }
  end
end

class BigWigComputationParallelizer < WigComputationParallelizer
  # Apparently, Ruby doesn't have super.method() calling...
  alias :super_run :run

  # Same as WigComputationParallelizer#run, but convert the final output to BigWig
  def run(wig, assembly, &block)
    super_run(wig, &block)
    
    # Convert the output Wig file to BigWig
    tmp_file = @output + '.tmp'
    Wig.to_bigwig(@output, tmp_file, assembly)
    
    # Delete the temporary intermediate Wig file by moving the BigWig on top of it
    FileUtils.move(tmp_file, @output)
  end
end

##
# Module for managing Wig processes in parallel
##
module WigParallelEnumerable
  @@pm = Parallel::ForkManager.new(2, {'tempdir' => Dir.tmpdir})
  @@default_chunk_size = 200_000
  
  # Set the total number of computation processes for all Wig files (default = 2)
  def self.max_threads=(n)
    @@pm = Parallel::ForkManager.new(n.to_i, {'tempdir' => Dir.tmpdir})
  end
  
  # Set the default chunk size
  def self.chunk_size=(n)
    @@default_chunk_size = n
  end
 
  # Run a given block for each chromosome
  # (parallel each)
  def p_each    
    self.chromosomes.each do |chr|
      @@pm.start(chr)
      yield(chr)
      @@pm.finish(0)
    end

    @@pm.wait_all_children
  end

  # Compute a given block for each chromosome
  # and return the results (parallel map)
  def p_map
    results = Array.new(self.chromosomes.length)
    
    @@pm.run_on_finish do |pid,exit_code,id,exit_signal,core_dump,data|
      results[id] = data['output']
    end

    self.chromosomes.each_with_index do |chr,i|
      @@pm.start(chr)
      result = yield(chr)
      @@pm.finish(0, {'output' => result})
    end
    
    @@pm.wait_all_children

    return results
  end
  
  # Iterate over all chromosomes in chunks
  def chunk_each(chunk_size = @@default_chunk_size)
    self.chromosomes.each do |chr|
      # Run in parallel processes managed by ForkManager
      @@pm.start(chr) and next 
      puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
      
      chunk_start = 1
      chr_length = self.chr_length(chr)
      while chunk_start < chr_length
        chunk_stop = [chunk_start+chunk_size-1, chr_length].min
        puts "Processing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
        yield(chr, chunk_start, chunk_stop)
        chunk_start = chunk_stop + 1
      end

      @@pm.finish(0)
    end

    @@pm.wait_all_children
  end
  
  # Compute a given block for all chromosomes in chunks
  # and inject the results (parallel inject)
  def chunk_map(initial_value, chunk_size = @@default_chunk_size)
    results = Array.new(self.chromosomes.length)
    
    @@pm.run_on_finish do |pid,exit_code,id,exit_signal,core_dump,data|
      results[id] = data['output']
    end

    self.chromosomes.each_with_index do |chr,i|
      # Run in parallel processes managed by ForkManager
      @@pm.start(i) and next   
      puts "\nProcessing chromosome #{chr}" if ENV['DEBUG']
      
      result = initial_value
      chunk_start = 1
      chr_length = self.chr_length(chr)
      while chunk_start < chr_length
        chunk_stop = [chunk_start+chunk_size-1, chr_length].min
        puts "Computing chunk #{chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
        
        result = yield(result, chr, chunk_start, chunk_stop)
        
        chunk_start = chunk_stop + 1
      end
      
      @@pm.finish(0, {'output' => result})
    end

    @@pm.wait_all_children
    
    return results
  end
end