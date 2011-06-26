require 'stats'

##
# Additional methods mixed in for EntryFile types that also
# have numeric data/values
##
module SpotFile
  # Get a Spot with a specific id
  def id(query_id)
    entries = Array.new
    skipped = 0
    File.grep(@data_file, query_id) do |line|
      begin
        entries << parse(line)
      rescue
        skipped += 1
      end
    end
    
    raise EntryFileError, "No spot with id #{query_id} in file #{File.basename(@data_file)}" if entries.length == 0
    raise EntryFileError, "More than one spot with id #{query_id} in file #{File.basename(@data_file)}" if entries.length > 1
    puts "Skipped #{skipped} invalid spots with id #{query_id}" if ENV['DEBUG']
    
    return entries.first
  end
  
  ##
  # STATISTICAL METHODS
  ##
  
  # The sum of the values of all spots
  def total
    # Cache for performance
    @total = self.map { |entry| entry.value }.sum if @total.nil?
    return @total
  end
  
  # The mean value of all spots
  def mean
    total.to_f / count
  end
  
  # The standard deviation of all spots
  def stdev(avg = self.mean)
    # Cache for performance
    if @stdev.nil?
      sum_of_deviances = self.map { |entry| (entry.value-avg)**2 }.sum
      @stdev = sum_of_deviances.to_f / count
    end
    
    return @stdev
  end
  
  ##
  # QUERY METHODS
  ##
  
  # Return a Contig of values for the given window
  def query(chr, start, stop)
    low = [start, stop].min
    high = [start, stop].max
    length = high - low + 1
    
    data = Contig.new(chr)
    
    self.each(chr, start, stop) do |spot|
      # Get the high and low spot coordinates, and clamp to the ends of the window
      low = [low, spot.low].max
      high = [spot.high, high].min
    
      for bp in low..high
        contig.set(bp, spot.value) unless spot.value.nil?
      end
    end
    
    return contig
  end
  
  ##
  # OUTPUT METHODS
  ##
  
  # Write this array to variableStep Wig format
  def to_wig(filename, assembly)    
    # Iterate over each chromosome, mapping all spots and averaging
    File.open(File.expand_path(filename), 'w') do |f|
      # TODO: should be rewritten to intelligently use step size
      f.puts UCSCTrackHeader.new(:type => 'wiggle_0').to_s
      
      self.chromosomes.each do |chr|
        # Skip if this chromosome is not in the specified assembly
        next unless assembly.include?(chr)
        
        # Allocate space for the new Wig chromosomes
        values = query(chr, 1, assembly[chr])
      
        # Write to output file
        f.puts values.to_s
      end
    end
  end
  
  # Write this array to BigWig format
  # By first writing to bedGraph, then calling BedGraph#to_bigwig
  def to_bigwig(filename, assembly)
    begin
      tmp_bedgraph = File.expand_path(filename + '.bedGraph')
      self.to_bedgraph(tmp_bedgraph)
    
      # bedGraph must be sorted to call bedGraphToBigWig
      tmp_sorted = tmp_bedgraph + '.sorted'
      File.sort(tmp_bedgraph, tmp_sorted, '-k1,1 -k2,2n')
      %x[ bedGraphToBigWig #{tmp_sorted} #{File.expand_path(assembly.len_file)} #{File.expand_path(filename)} ]
    rescue
      raise "Error converting Array to BigWig!"
    ensure
      # Delete the temporary intermediate files
      File.delete(tmp_bedgraph) if File.exist?(tmp_bedgraph)
      File.delete(tmp_sorted) if File.exist?(tmp_sorted)
    end
  end
end