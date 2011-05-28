# Use GSL-based statistical methods for now
# Other implementations include:
# - RStats: relying on rs-ruby gem / R
# - NativeStats: all methods implemented in Ruby code
# - GSLStats: relying on rb-gsl / GSL (GNU Scientific Library)
# To employ these alternatives, change the include

require 'stats/gsl_stats'
require 'stats/native_stats'

module Enumerable
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

class Array
    include NativeStats
end