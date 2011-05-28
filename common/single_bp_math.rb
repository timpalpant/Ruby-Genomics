#
#  SingleBPMath.rb
#  BioRuby
#
#  Some operations for genomic math
#
#  Created by Timothy Palpant on 3/28/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'stats'

module SingleBPMath
	# Return the total number of values in this GenomicData
  def num_values
  	self.collect { |chr_id,data| data.length }.sum.to_i
  end
	
	# Return the sum of all values
  def total
    self.collect { |chr_id,values| values.sum }.sum
  end

  # Return the mean of all values
  def mean
  	total.to_f / num_values
  end
	
	# Return the standard deviation of all values
  def stdev(mean = self.mean, num_values = self.num_values)
  	tss = self.map { |chr_id,values| values.tss(mean) }.sum
  	Math.sqrt(tss / num_values)
  end
end