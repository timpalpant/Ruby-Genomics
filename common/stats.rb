# Use native implementations of statistical methods for now
# Other implementations include:
# - RStats: relying on rs-ruby gem / R
# - NativeStats: all methods implemented in Ruby code
# - GSLStats: relying on rb-gsl / GSL (GNU Scientific Library)

require 'stats/native_stats'

class Array
  include NativeStats

	# Return true if any elements are true
	def any?
		self.each { |e| return true if e }
		return false
	end
	
	# Return true if all elements are true
	def all?
		self.each { |e| return false unless e }
		return true
	end
end