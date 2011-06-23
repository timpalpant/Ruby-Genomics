require 'stats'

##
# Methods for basic statistics of SpotArrays / SpotFiles
##
module SpotArrayMath
  # The total number of spots in this SpotArray
  def num_spots
    self.chromosomes.map { |chr| self[chr].length }.sum
  end
  
  # The sum of the values of all spots
  def total
    self.chromosomes.map { |chr| self[chr].collect { |spot| spot.value }.sum }.sum
  end
  
  # The mean value of all spots
  def mean
    total.to_f / num_spots
  end
  
  # The standard deviation of all spots
  def stdev(mean = self.mean)
    self.chromosomes.map { |chr| self[chr].map { |spot| (spot.value-mean)**2 }.sum }.sum / num_spots.to_f
  end
end