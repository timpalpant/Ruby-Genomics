require 'stats'

##
# Sugar for computing basic statistics on Wig files
##
module WigMath
  # Number of values in the Wig file
  def num_values
    self.chromosomes.map { |chr| chr_length(chr) }.sum
  end
  
  # The sum of all values
  def total
    chr_totals = self.chunk_map(0) do |sum,chr,start,stop|
      sum + query(chr,start,stop).sum
    end

    return chr_totals.sum
  end
  
  # The mean of all values
  def mean
    total / num_values
  end
  
  # The standard deviation of all values
  def stdev(avg = self.mean)
    chr_deviances = self.chunk_map(0) do |sum,chr,start,stop|
      sum + query(chr,start,stop).map { |elem| (elem-avg)**2 }.sum
    end
    
    Math.sqrt(chr_deviances.sum / num_values)
  end
end