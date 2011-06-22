#
#  fft.rb
#  BioRuby
#  Compute the DFT and related spectra of an Array with GSL
#
#  Created by Timothy Palpant on 4/8/11.
#  Copyright 2011 UNC. All rights reserved.
#

# From the GSL documentation:
# The functions use the FFTPACK storage convention for half-complex sequences. In this convention the half-complex transform of a real sequence is stored with frequencies in increasing order, starting at zero, with the real and imaginary parts of each frequency in neighboring locations. When a value is known to be real the imaginary part is not stored. The imaginary part of the zero-frequency component is never stored. It is known to be zero (since the zero frequency component is simply the sum of the input data (all real)). For a sequence of even length the imaginary part of the frequency n/2 is not stored either, since the symmetry z_k = z_{N-k}^* implies that this is purely real too.

# The storage scheme is best shown by some examples. The table below shows the output for an odd-length sequence, n=5. The two columns give the correspondence between the 5 values in the half-complex sequence returned by gsl_fft_real_transform, halfcomplex[] and the values complex[] that would be returned if the same real input sequence were passed to gsl_fft_complex_backward as a complex sequence (with imaginary parts set to 0),

# complex[0].real  =  halfcomplex[0] 
# complex[0].imag  =  0
# complex[1].real  =  halfcomplex[1] 
# complex[1].imag  =  halfcomplex[2]
# complex[2].real  =  halfcomplex[3]
# complex[2].imag  =  halfcomplex[4]
# complex[3].real  =  halfcomplex[3]
# complex[3].imag  = -halfcomplex[4]
# complex[4].real  =  halfcomplex[1]
# complex[4].imag  = -halfcomplex[2]
# The upper elements of the complex array, complex[3] and complex[4] are filled in using the symmetry condition. The imaginary part of the zero-frequency term complex[0].imag is known to be zero by the symmetry.

# The next table shows the output for an even-length sequence, n=5 In the even case both the there are two values which are purely real,

# complex[0].real  =  halfcomplex[0]
# complex[0].imag  =  0
# complex[1].real  =  halfcomplex[1] 
# complex[1].imag  =  halfcomplex[2] 
# complex[2].real  =  halfcomplex[3] 
# complex[2].imag  =  halfcomplex[4] 
# complex[3].real  =  halfcomplex[5] 
# complex[3].imag  =  0 
# complex[4].real  =  halfcomplex[3] 
# complex[4].imag  = -halfcomplex[4]
# complex[5].real  =  halfcomplex[1] 
# complex[5].imag  = -halfcomplex[2] 
# The upper elements of the complex array, complex[4] and complex[5] are filled in using the symmetry condition. Both complex[0].imag and complex[3].imag are known to be zero.

require 'gsl'
include GSL

module FFT
  
  # Compute the power spectrum
  def power_spectrum
    # Apply the DFT to the original input data
    # and compute the power spectrum from the returned complex Fourier coefficients
    # Subset the returned results to strip off the first (0 frequency) and last (N/2 frequency)
    # because of how they are returned (see above)
    return self.to_gslv.fft.subvector(1, self.length-2).to_complex2.abs2.to_a
  end
  
  # Compute a spectrogram with window size w
  # TODO: Ends?
  def spectrogram(w)
    v = self.to_gslv
    half_window = w/2
    output = Matrix[half_window, self.length]
    
    for i in half_window...self.length-half_window
      low = i - half_window
      high = i + half_window

      output.set_col(i, v[low..high].power_spectrum)
    end
    
    return output
  end
  
end

class Array
  include FFT
end

class Vector
  include FFT

  def to_gslv
    self
  end
end
