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
		self.sum.to_f / numel unless numel.zero?
	end
	
	def variance
		avg = self.mean
		sum_of_deviance = self.compact.inject(0) { |sum, elem| sum + (elem-avg)**2 }
		return sum_of_deviance / self.compact.size
	end
	
  def stdev
  	Math.sqrt(variance)
  end	
  
	# Struct to contain index-value pairs for index_sort
	IndexValuePair = Struct.new(:index, :value)
	
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
	def index_sort
		sorted_indices  = Array.new(self.length)
		self.each_with_index do |elem,i|
			sorted_indices[i] = IndexValuePair.new(i,elem)
		end
		
		# Sort by element but keep indices
		sorted_indices.sort { |e1,e2| e1.value <=> e2.value }.map { |elem| elem.index }
	end

	def zscore
		avg = self.mean     
		standard_dev = self.stdev
		self.map { |elem| (elem-avg)/standard_dev unless elem.nil? }
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

	# Smooth the array with a window (mean)
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
end