require 'bio/wig'

module Bio
  class WigFile
    # Compute some transformation on the values in this file, and output the result to disk
    # Should take a block that accepts coordinates: chr, start, stop
    # and returns a Contig of values for those coordinates
    # This transformed Contig will then be written to disk, and all of the transformed chunks
    # will be concatenated into a single, transformed output (text) Wig file
    def transform(output_file, assembly = nil, opts = {})
      # Write the output file header
      header_file = output_file+'.header'
      File.open(header_file, 'w') do |f|
        f.puts @track_header.to_s
      end

      # Keep track of all the temporary intermediate files (header first)
      tmp_files = [header_file]
      @contigs_index.each { |contig_info| tmp_files << "#{output_file}.#{contig_info.chr}.#{contig_info.start}" }
  
      # Iterate by chromosome
      @contigs_index.p_each(opts) do |contig_info|
        puts "\nProcessing contig #{contig_info}" if ENV['DEBUG']

        # Write the header
        chr_temp_file = "#{output_file}.#{contig_info.chr}.#{contig_info.start}"
        File.open(chr_temp_file, 'w') do |f|
          f.puts "variableStep chrom=#{contig_info.chr} span=1"
        end
      
        chunk_start = contig_info.start
        while chunk_start <= contig_info.stop
          chunk_stop = [chunk_start+CHUNK_SIZE-1, contig_info.stop].min
          puts "Processing chunk #{contig_info.chr}:#{chunk_start}-#{chunk_stop}" if ENV['DEBUG']
        
          output = yield(contig_info.chr, chunk_start, chunk_stop)
        
          # Write this chunk to disk
          File.open(chr_temp_file, 'a') do |f|
            f.puts output.to_variable_step
          end
        
          chunk_start = chunk_stop + 1
        end
      end
      
      # Concatenate all of the temp file pieces into the final output
      File.cat(tmp_files, output_file)
    ensure
      # Delete the individual temp files created by each process
      tmp_files.each { |filename| File.delete(filename) if File.exist?(filename) }
    end
  end
end