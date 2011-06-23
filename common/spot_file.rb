require 'stats'

##
# Additional methods mixed in for EntryFile types that also
# have numeric data/values
##
module SpotFile
  ##
  # STATISTICAL METHODS
  ##
  
  # The sum of the values of all spots
  def total
    # Cache for performance
    if @total.nil?
      @total = 0
      self.each { |entry| @total += entry.value }
    end
    
    return @total
  end
  
  # The mean value of all spots
  def mean
    total.to_f / count
  end
  
  # The standard deviation of all spots
  def stdev(mean = self.mean)
    # Cache for performance
    if @stdev.nil?
      sum_of_deviances = 0
      self.each { |entry| sum_of_deviances = (entry.value-mean)**2 }
      @stdev = sum_of_deviances / count.to_f
    end
    
    return @stdev
  end
  
  ##
  # QUERY METHODS
  ##
  
  # Return values for the given window, with single-bp resolution (even if inferred)
  def query(chr, start, stop)
    low = [start, stop].min
    high = [start, stop].max
    length = high - low + 1
    
    total = Array.new(length, 0)
    count = Array.new(length, 0)
    
    self.each(chr, start, stop) do |spot|
      # Get the high and low spot coordinates, and clamp to the ends of the window
      low = [low, spot.low].max
      high = [spot.high, high].min
    
      for bp in spot.low..spot.high
        total[bp-low] += spot.value unless spot.value.nil?
        count[bp-low] += 1
      end
    end
    
    # Map base pairs without data to nil, and take the mean of overlapping probes
    avg = Array.new(length) do |i|
      if count[i] > 0
        total[i].to_f / count[i]
      else
        nil
      end
    end
    
    # Allow Crick querying
    avg.reverse! if start > stop
    return avg.to_contig(chr, start, 1, 1)
  end
  
  ##
  # OUTPUT METHODS
  ##
  
  # Write this array to Wig format
  # Construct a Wig that is the data from all spots
  # averaged, if multiple spots cover a given base
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
  def to_bigwig(filename, assembly)
    self.to_wig(filename, assembly)
    
    tmp_file = filename + '.tmp'
    WigFile.to_bigwig(filename, tmp_file, assembly)
    FileUtils.move(tmp_file, filename)
  end
end