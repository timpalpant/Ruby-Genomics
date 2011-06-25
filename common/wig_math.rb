require 'stats'

##
# Sugar for computing basic statistics on Wig files
##
module WigMath
  # Number of values in the Wig file
  def num_values
    @contigs_index.map { |chr,start,length| length }.sum
  end
  
  # The sum of all values
  def total
    sum = 0
    self.each_chunk do |chunk|
      sum += chunk.sum
    end

    return sum
  end
  
  # The mean of all values
  def mean
    total.to_f / num_values
  end
  
  # The standard deviation of all values
  def stdev(avg = self.mean)
    deviances = 0
    self.each_chunk do |chunk|
      deviances += chunk.map { |elem| (elem-avg)**2 }.sum
    end
    
    Math.sqrt(chr_deviances.sum / num_values)
  end
end