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
    raise 'Cannot convolve vectors of unequal length!' if a.length != b.length
  
    # Convolution is multiplication in frequency space
    v1 = a.to_gslv.fft.subvector(1, self.length-2).to_complex2
    v2 = b.to_gslv.fft.subvector(1, self.length-2).to_complex2
    
		return (v1 * v2).ifft.to_a
	end
  
  # Should only be used if computing a gaussian a few times,
  # otherwise, the filter should be Fourier transformed once
  # and then it can be reused
  def self.gaussian_filter(a, sdev, window_m = 3)
    g = Filter.gaussian(sdev, window_m)
    
    # Pad the gaussian vector with zeros to make it the same size as a
    padded = Vector[a.length]
    padded[a.length/2-g.length/2..a.length/2+g.length/2] = gaussian
    
    # Convolve the array with the Gaussian filter
    convolve(a, padded)
	end
end

class Filter < Array
  # Return a Gaussian vector (of length 2*sdev*window_m)
  def self.gaussian(sdev, window_m = 3)
    half_window = sdev*window_m
    gaussian = Vector[2*half_window]
    #coeff = 1 / (sdev*Math.sqrt(2*Math::PI))
    for x in -half_window..half_window
      gaussian[x+half_window] = Math.exp(-((x**2)/(2*(sdev**2))))
    end
    gaussian /= gaussian.sum
    
    return self.new(gaussian.to_a)
  end
  
  # Compute the Fourier transform of this filter (returns a GSL::ComplexVector)
  def fft
    self.to_gslv.fft.subvector(1, self.length-2).to_complex2
  end
end