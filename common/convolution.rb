#
#  convolution.rb
#  BioRuby
#  Compute convolutions of Arrays with GSL
#
#  Created by Timothy Palpant on 6/1/11.
#  Copyright 2011 UNC. All rights reserved.
#

# Compute efficient convolutions with the FFT

require 'gsl'
include GSL

module Math
	# Convolve Arrays a and b
	def self.convolve(a,b)
    raise 'Cannot convolve Arrays of unequal length!' if a.length != b.length
  
    # Convolution is multiplication in frequency space
    v1 = a.to_gslv.fft
    v2 = b.to_gslv.fft
 
		return (v1 * v2).ifft.to_complex.fftshift.real.to_a
	end
  
  # Should only be used if computing a gaussian a few times,
  # otherwise, the filter should be Fourier transformed once
  # and then it can be reused
  def self.gaussian_filter(a, sdev, window_m = 3)
    g = Filter.gaussian(sdev, window_m)
    
    # Pad the gaussian vector with zeros to make it the same size as a
    padded = Vector[a.length]
    padded[a.length/2-g.length/2..a.length/2+g.length/2] = g
    
    # Convolve the array with the Gaussian filter
    convolve(a, padded)
	end
end

class Filter < Array
  # Return a Gaussian vector (of length 2*sdev*window_m)
  def self.gaussian(sdev, window_m = 3)
    half_window = sdev*window_m
    g = Vector[2*half_window+1]
    #coeff = 1 / (sdev*Math.sqrt(2*Math::PI))
    for x in -half_window..half_window
      g[x+half_window] = Math.exp(-((x**2)/(2*(sdev**2))))
    end
    g /= g.sum
    
    return self.new(g.to_a)
  end
end

class Vector
  def to_gslv
    self
  end
end
