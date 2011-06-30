# Use native implementations of statistical methods for now
# Other implementations include:
# - RStats: relying on rs-ruby gem / R
# - NativeStats: all methods implemented in Ruby code
# - GSLStats: relying on rb-gsl / GSL (GNU Scientific Library)

require 'stats/native_stats'

class Array
  include NativeStats
end
