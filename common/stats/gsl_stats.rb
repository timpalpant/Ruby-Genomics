require 'gsl'
include GSL

##
# Add methods to compute descriptive statistics on GSL::Vectors
##

class Vector
  def gaussian_smooth(sdev,window_size)
    half_window = window_size*sdev
    
    # Generate the gaussian vector (of length window_size)
    gaussian = Vector[2*half_window+1]
    coeff = 1 / (sdev*Math.sqrt(2*Math::PI))
    for x in -half_window..half_window
      gaussian[x+half_window] = coeff * Math.exp(-((x**2)/(2*(sdev**2))))
    end
    
    smooth = Vector[self.length]
    for i in half_window...self.length-half_window
      window = self[i-half_window..i+half_window]
      smooth[i] = (window * gaussian).sum
			# Hack to force GC
			GC.start if i % 500 == 0 # be a little less aggressive about GC
    end
  
    return smooth
  end
end