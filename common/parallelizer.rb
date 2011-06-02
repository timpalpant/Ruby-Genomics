#
#  parallelizer.rb
#  BioRuby
#  Parallelize computations across chromosomes with ForkManager
#
#  Created by Timothy Palpant on 6/1/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'forkmanager'

class Parallelizer
  def initialize(max_threads = 2)
    @pm = Parallel::ForkManager.new(max_threads)
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
      f.puts Wig.track_header(@output, @output)
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
        f.puts Wig.fixed_step(chr) + ' start=1 step=1 span=1'
      end
      
      chunk_start = 1
      chr_length = wig.chr_length(chr)
      while chunk_start < chr_length
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
    
    # Return a handle to the output file
    #WigFile.new(@output)
  end
end