# Loads and models a Position Specific Scoring Matrix or Position Weight Matrix (PWM)
# Also contains methods to score a sequence of the appropriate (PWM) length,
# i.e. compute the binding affinity energy
class PWM < Hash
  
  # Header summary, minimum possible energy (corresponding to the consensus sequence), 
  # and maximum possible energy
  attr_accessor :summary, :minE, :maxE
  
  def self.load(pwm_file, genome)
    # Store the background GC-content across the entire genome
    background = {'a' => genome.a_content.to_f, 'c' => genome.c_content.to_f, 'g' => genome.g_content.to_f, 't' => genome.t_content.to_f }
    
    # Initialize the hash with a default value of 0
    pwm = self.new(0)
    
    File.open(pwm_file) do |f|
      # Load summary line
      pwm.summary = f.gets
      
      # Read in frequencies and map to log(freq/background) for energy scoring
      # TODO: Ensure normalization for non-standard PWMs
      # A frequencies
      pwm['a'] = f.gets.chomp.split("\t")[1..-1].map { |value| Math.log(value.to_f / background['a']) }
      # C frequencies
      pwm['c'] = f.gets.chomp.split("\t")[1..-1].map { |value| Math.log(value.to_f / background['c']) }
      # G frequencies
      pwm['g'] = f.gets.chomp.split("\t")[1..-1].map { |value| Math.log(value.to_f / background['g']) }
      # T frequencies
      pwm['t'] = f.gets.chomp.split("\t")[1..-1].map { |value| Math.log(value.to_f / background['t']) }
    end
    
    # BPs in the PWM should be the same length
    raise "Invalid PWM!" if pwm['a'].length != pwm['c'].length or pwm['c'].length != pwm['g'].length or pwm['g'].length != pwm['t'].length
    
    # Store minimum possible energy for this PWM (the consensus sequence)
    max_score = 0
    for i in 0...pwm.length
      max_score += pwm.at(i).values.max
    end
    pwm.minE = -max_score
    
    # Store maximum possible energy for this PWM (the least likely to bind)
    min_score = 0
    for i in 0...pwm.length
      min_score += pwm.at(i).values.min
    end
    pwm.maxE = -min_score
      
    return pwm
  end
  
  def to_s
    @summary
  end
  
  # Return a hash of the values at a given index
  def at(index)
    raise "Index out of range!" if index < 0 or index >= self.length
    
    bp_weights = Hash.new
    self.each_key { |bp| bp_weights[bp] = self[bp][index] }
    return bp_weights
  end
  
  # PWM shouldn't have been created with an irregular matrix, so we'll assume
  # the 'a' length is the sequence length
  # TODO: Should be a little more thorough
  def length
    self['a'].length
  end
  
  # Compute an energy score for the given sequence
  def score(sequence)
    raise "Sequence length #{sequence.length} does not match PWM length #{self.length}!" if sequence.length != self.length
    
    pwm_score = 0
    sequence.downcase.split(//).each_with_index do |bp,i|
      pwm_score += self[bp][i]
    end
    
    return -pwm_score - @minE
  end
  
end