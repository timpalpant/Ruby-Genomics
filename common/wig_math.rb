require 'stats'

##
# Sugar for computing basic statistics on Wig files
##
module WigMath
  # Number of values in the Wig file
  def num_bases
    # Cache for performance
    compute_stats if @num_bases.nil?
    return @num_bases
  end
  
  # The sum of all values
  def total
    # Cache for performance
    compute_stats if @total.nil?
    return @total
  end
  
  # The mean of all values
  def mean  
    # Cache for performance
    compute_stats if @mean.nil?
    return @mean unless num_bases == 0
  end
  
  # The standard deviation of all values
  def stdev(avg = self.mean)
    compute_stats if @stdev.nil?
    return @stdev unless num_bases == 0
  end
  
  private
  
  # Compute the coverage, total, and stdev in a single iteration
  def compute_stats
    @num_bases = 0
    @total = 0
    sum_of_squares = 0.0
    
    self.each_chunk do |chunk|
      @num_bases += chunk.coverage
      @total += chunk.sum
      sum_of_squares += chunk.values.compact.map { |elem| elem**2 }.sum
    end
    
    @mean = @total.to_f / @num_bases
    variance = (sum_of_squares - @total*@mean) / @num_bases
    @stdev = Math.sqrt(variance)
  end
end
