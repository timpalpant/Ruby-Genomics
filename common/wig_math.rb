require 'stats'

##
# Sugar for computing basic statistics on Wig files
##
module WigMath
  # Number of values in the Wig file
  def num_bases
    count = 0
    self.each_chunk do |chunk|
      count += chunk.coverage
    end
    
    return count
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
    total.to_f / num_bases
  end
  
  # The standard deviation of all values
  def stdev(avg = self.mean)
    deviances = 0.0
    self.each_chunk do |chunk|
      deviances += chunk.to_a.compact.map { |elem| (elem-avg)**2 }.sum
    end
    
    Math.sqrt(deviances / num_bases)
  end
end