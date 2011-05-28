##
# Static math functions that are either missing or inefficient in Ruby
##
module Math
  def self.max(a, b)
		a > b ? a : b unless a.nil? or b.nil?
	end

	def self.min(a, b)
		a < b ? a : b unless a.nil? or b.nil?
	end
end