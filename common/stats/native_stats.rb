##
# Methods to compute descriptive statistics
# using native Ruby implementations
##
module NativeStats
  def sum
    self.compact.inject { |sum, elem| sum + elem }
  end
  
  def mean
    numel = self.compact.length
    self.sum.to_f / numel unless numel == 0
  end
  
  def variance(avg = self.mean)
    compacted = self.compact
    return nil if compacted.length == 0
    sum_of_deviance = compacted.map { |elem| (elem-avg)**2 }.sum
    return sum_of_deviance / compacted.length
  end
  
  def stdev(avg = self.mean)
    return nil if self.compact.length == 0
    Math.sqrt(variance(avg))
  end 
  
  def median
    sorted = self.compact.sort
    if sorted.length == 0
      nil
    elsif sorted.length.odd?
      # Median is the middle value
      sorted[sorted.length/2]
    else
      # Median is the average of the middle two values
      (sorted[sorted.length/2-1] + sorted[sorted.length/2]) / 2.0
    end
  end
  
  # Lower quartile is the median value of the elements less than the median
  def lower_quartile
    sorted = self.sort
    midpoint = sorted.length/2 - 1 
    sorted[0..midpoint].median
  end
  
  # Upper quartile is the median value of the elements greater than the median
  def upper_quartile
    sorted = self.sort
    midpoint = if sorted.length.even?
      sorted.length/2
    else
      sorted.length/2 + 1
    end
    sorted[midpoint..-1].median
  end

  # Sort an array and return the index (like in Matlab)
  def sort_index
    indices = (0...self.length).to_a
    
    # Sort the indices by their element (but keep indices)
    indices.sort { |e1,e2| self[e1] <=> self[e2] }
  end

  def zscore
    avg = self.mean     
    sdev = self.stdev
    self.map { |elem| (elem-avg)/sdev unless elem.nil? }
  end

  # Returns a hash of objects and their frequencies within array.
  def freq                                 
    k = Hash.new(0)
    self.each { |x| k[x] += 1 }
    return k
  end

  # Given two arrays a and b, a^b returns a new array of objects *not* found in the union of both.
  def ^(other)                             
  (self | other) - (self & other)
  end

  # Return the value of the pth percentile
  def percentile(p)
    self.sort[(p * self.length).ceil - 1]
  end

  # Returns the frequency of x within array.
  def freq(x)                              
    count[x]
  end

  # Returns highest count of any object within array.
  def maxcount                              
    count.values.max
  end

  # Returns lowest count of any object within array.
  def mincount                              
    count.values.min
  end

  # Returns a new array of object(s) with x highest count(s) within array.
  def outliers(x)                           
    h = count                                                              
    min = count.values.uniq.sort.reverse.first(x).min
    h.delete_if { |x,y| y < min }.keys.sort
  end

  # Smooth the array with a moving average window (mean)
  def window_smooth(window_size)
    return self if window_size == 1
    
    half_window = window_size / 2   
    moving_sum = self[0...half_window].sum.to_f
    nil_values = self[0...half_window].map { |elem| elem.nil? ? 1 : 0 }.sum
    
    moving_average = Array.new
    self.each_index do |i|
      lose = [i-half_window-1, -1].max
      gain = [self.length, i+half_window].min
      
      if self[lose].nil?
        nil_values -= 1
      elsif lose != -1
        moving_sum -= self[lose]
      end
      
      if self[gain].nil?
        nil_values += 1
      elsif gain != self.length
        moving_sum += self[gain]
      end
      
      avg = moving_sum / (gain-lose-nil_values)
      # Replace NaN and Infinity with nil
      moving_average << (avg.finite? ? avg : nil)
    end
    
    return moving_average
  end
  
  # Smooth the Array with a Gaussian filter
  def gaussian_smooth(sdev, window_size)
    half_window = window_size*sdev
    
    # Generate the gaussian vector (of length window_size)
    gaussian = Array.new(2*half_window+1, 0)
    coeff = 1 / (sdev*Math.sqrt(2*Math::PI))
    for x in -half_window..half_window
      gaussian[x+half_window] = coeff * Math.exp(-((x**2)/(2*(sdev**2))))
    end
    
    smooth = Array.new(self.length, 0)
    i = half_window
    self.each_cons(2*half_window+1) do |window|
      for j in 0...window.length
        smooth[i] += window[j] * gaussian[j]
      end

      i += 1
    end

    return smooth
  end
end
