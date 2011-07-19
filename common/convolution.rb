#
#  convolution.rb
#  BioRuby
#  Compute convolutions with fftw3
#  Good resource: http://morse.cs.byu.edu/450/lectures/lect21/fft.slides.printing.6.pdf
#
#  Created by Timothy Palpant on 6/1/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'fftw3'

module Math
  # Convolve the real NArrays a and b
  def self.convolve(a,b)
    raise "Cannot convolve NArrays of unequal length!" if a.length != b.length
  
    # Convolution is multiplication in frequency space:
    # Take the DFT
    f1 = FFTW3.fft(a)
    f2 = FFTW3.fft(b)
    # Multiply in frequency space, and invert the DFT
    FFTW3.ifft(f1 * f2).real / a.length
  end
  
  # Apply the filter f to the NArray a
  # Convenience method for applying filters without having to figure out padding
  def self.filter(a, f)
    convolve(a, f.fft(a.length))[0...a.length]
  end
end

class Filter < NArray
  # Return a box filter (of length n)
  def self.box(n)
    self.float(n).fill(1.0 / n)
  end
  
  # Return a Gaussian filter (of length 2*sdev*window_m)
  def self.gaussian(sdev, window_m = 3)
    half_window = sdev*window_m
    g = self.float(2*half_window+1)
    for x in -half_window..half_window
      g[x+half_window] = Math.exp(-((x**2)/(2*(sdev**2))))
    end
    g /= g.sum
    
    return g
  end
  
  def fft(n)
  
  end
end
