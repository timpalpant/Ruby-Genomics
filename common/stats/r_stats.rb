require 'rsruby'

##
# Expose R functionality in a singleton module
##
module R
  
  @@R = RSRuby.instance

  LSFit = Struct.new(:slope, :intercept)
  
  def self.assign(name, data)
    @@R.assign(name, data)
  end
  
  def self.eval(statement)
    @@R.eval_R(statement)
  end
  
  def self.lsfit(x, y)
    # TODO: Find a cleaner way to accomplish curve fitting
    @@R.assign('x', x)
    @@R.assign('y', y)
    fit = @@R.lm('x ~ y')
    @@R.remove('x', 'y')
    
    return LSFit.new(fit['coefficients']['y'], fit['coefficients']['(Intercept)'])
  end
  
  def self.correlation(x, y)
    @@R.cor(x, y)
  end
end


##
# Add methods to the Array class from R
##
module RStats
  @@R = RSRuby.instance
  
  def mean
    @@R.mean(self)
  end
  
  # Returns a moving average with window_size *around* each element
  def smooth(window_size)
    return self if window_size == 1
    
    half_window = window_size / 2
    
    moving_average = Array.new
    self.each_index do |i|
      start = [i-half_window, 0].max
      stop = [self.length-1, i+half_window].min
      
      moving_average << self[start..stop].mean
    end
    
    return moving_average
  end
end