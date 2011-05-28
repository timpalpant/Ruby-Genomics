##
# Add methods to Float and Fixnum classes to output a fixed number
# of significant digits by using an optional argument
##
class Float
  alias_method :orig_to_s, :to_s
  def to_s(arg = nil)
    if arg.nil?
      orig_to_s
    else
      # Outputs a fixed number of decimal places
      #sprintf("%.#{arg}f", self)
      # Outputs a fixed number of digits (significant figures)
      "%.#{arg}g" % self
    end
  end
end

class Fixnum
  # Added so that to_s can be called on Floats or Fixnums
  # without being interpreted as a radix (see RubyDoc for Integer#to_s(radix))
  alias_method :orig_to_s, :to_s
  def to_s(arg = nil)
    orig_to_s()
  end
end