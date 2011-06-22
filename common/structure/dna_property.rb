require 'bio'

require 'orchid'
require 'pwm'
require 'stats'

# soap4r is broken in Ruby 1.9.2
# TODO: Implement GenomeAtlasAPI with Savon
#require 'structure/genome_atlas_api'


##
# Compute DNA structural properties on a DNA Sequence (add methods to BioRuby)
##
class Bio::Sequence::NA

  # Singleton representing the GenomeAtlas SOAP web service
  # soap4r is broken in Ruby 1.9.2
  #@@service = GenomeAtlasAPI.new
  
  # Dinucleotide structural constants used in local calculations (from http://srs6.bionet.nsc.ru/srs6bin/cgi-bin/wgetz?-page+FieldInfo+-id+6F1Iv1SJDVl+-lib+PROPERTY+-bf+PropertyName
  # Originally from Orson et al. 1996 (?)
  @@PROPELLER = {'aa'=>-17.3, 'ac'=>-6.7, 'ag'=>-14.3, 'at'=>-16.9, 'ca'=>-8.6, 'cc'=>-12.8, 'cg'=>-11.2, 'ct'=>-14.3, 'ga'=>-15.1, 'gc'=>-11.7, 'gg'=>-12.8, 'gt'=>-6.7, 'ta'=>-11.1, 'tc'=>-15.1, 'tg'=>-8.6, 'tt'=>-17.3 }
  @@SLIDE = {'aa'=>-0.03, 'ac'=>-0.13, 'ag'=>0.47, 'at'=>-0.37, 'ca'=>1.46, 'cc'=>0.6, 'cg'=>0.63, 'ct'=>0.47, 'ga'=>-0.07, 'gc'=>0.29, 'gg'=>0.6, 'gt'=>-0.13, 'ta'=>0.74, 'tc'=>-0.07, 'tg'=>1.46, 'tt'=>-0.03 }
  @@STACKING_ENERGY = {'aa'=>-5.37, 'ac'=>-10.51, 'ag'=>-6.78, 'at'=>-6.57, 'ca'=>-6.57, 'cc'=>-8.26, 'cg'=>-9.61, 'ct'=>-6.78, 'ga'=>-9.81, 'gc'=>-14.59, 'gg'=>-8.26, 'gt'=>-10.51, 'ta'=>-3.82, 'tc'=>-9.81, 'tg'=>-6.57, 'tt'=>-5.37 }
  @@FREE_ENERGY = {'aa'=>-1.2, 'ac'=>-1.5, 'ag'=>-1.5, 'at'=>-0.9, 'ca'=>-1.7, 'cc'=>-2.1, 'cg'=>-2.8, 'ct'=>-1.5, 'ga'=>-1.5, 'gc'=>-2.3, 'gg'=>-2.1, 'gt'=>-1.5, 'ta'=>-0.9, 'tc'=>-1.5, 'tg'=>-1.7, 'tt'=>-1.2 }
  @@ENTHALPY = {'aa'=>-8.0, 'ac'=>-9.4, 'ag'=>-6.6, 'at'=>-5.6, 'ca'=>-8.2, 'cc'=>-10.9, 'cg'=>-11.8, 'ct'=>-6.6, 'ga'=>-8.8, 'gc'=>-10.5, 'gg'=>-10.9, 'gt'=>-9.4, 'ta'=>-6.6, 'tc'=>-8.8, 'tg'=>-8.2, 'tt'=>-8.0 }
  @@ENTROPY = {'aa'=>-21.9, 'ac'=>-25.5, 'ag'=>-16.4, 'at'=>-15.2, 'ca'=>-21.0, 'cc'=>-28.4, 'cg'=>-29.0, 'ct'=>-16.4, 'ga'=>-23.5, 'gc'=>-26.4, 'gg'=>-28.4, 'gt'=>-25.5, 'ta'=>-18.4, 'tc'=>-23.5, 'tg'=>-21.0, 'tt'=>-21.9 }
  @@MELTING = {'aa'=>54.5, 'ac'=>97.73, 'ag'=>58.42, 'at'=>57.02, 'ca'=>54.71, 'cc'=>85.97, 'cg'=>72.55, 'ct'=>58.42, 'ga'=>86.44, 'gc'=>136.12, 'gg'=>85.97, 'gt'=>97.73, 'ta'=>36.73, 'tc'=>86.44, 'tg'=>54.71, 'tt'=>54.5 }
  @@CLASH_STRENGTH = {'aa'=>0.64, 'ac'=>0.95, 'ag'=>2.53, 'at'=>1.68, 'ca'=>0.8, 'cc'=>1.78, 'cg'=>2.42, 'ct'=>2.53, 'ga'=>0.03, 'gc'=>0.22, 'gg'=>1.78, 'gt'=>0.95, 'ta'=>0, 'tc'=>0.03, 'tg'=>0.8, 'tt'=>0.64}
  @@MAJOR_GROOVE_MOBILITY = {'aa'=>1.18, 'ac'=>1.06, 'ag'=>1.06, 'at'=>1.12, 'ca'=>1.06, 'cc'=>0.99, 'cg'=>1.02, 'ct'=>1.04, 'ga'=>1.08, 'gc'=>0.98, 'gg'=>1, 'gt'=>1.02, 'ta'=>1.07, 'tc'=>1.03, 'tg'=>1.03, 'tt'=>1.09}
  @@MAJOR_GROOVE_SIZE = {'aa'=>3.98, 'ac'=>3.98, 'ag'=>4.7, 'at'=>4.7, 'ca'=>3.98, 'cc'=>3.98, 'cg'=>4.7, 'ct'=>4.7, 'ga'=>3.26, 'gc'=>3.26, 'gg'=>3.98, 'gt'=>3.98, 'ta'=>3.26, 'tc'=>3.26, 'tg'=>3.98, 'tt'=>3.98}
  @@MINOR_GROOVE_MOBILITY = {'aa'=>1.04, 'ac'=>1.1, 'ag'=>1.09, 'at'=>1.02, 'ca'=>1.16, 'cc'=>1.27, 'cg'=>1.25, 'ct'=>1.16, 'ga'=>1.12, 'gc'=>1.17, 'gg'=>1.25, 'gt'=>1.11, 'ta'=>1.05, 'tc'=>1.2, 'tg'=>1.23, 'tt'=>1.04}
  @@MINOR_GROOVE_SIZE = {'aa'=>2.98, 'ac'=>3.26, 'ag'=>3.98, 'at'=>3.26, 'ca'=>3.7, 'cc'=>3.98, 'cg'=>4.7, 'ct'=>3.98, 'ga'=>2.98, 'gc'=>3.26, 'gg'=>3.98, 'gt'=>3.26, 'ta'=>2.7, 'tc'=>2.98, 'tg'=>3.7, 'tt'=>2.98}
  @@RISE = {'aa'=>3.16, 'ac'=>3.41, 'ag'=>3.63, 'at'=>3.89, 'ca'=>3.23, 'cc'=>4.08, 'cg'=>3.6, 'ct'=>3.63, 'ga'=>3.47, 'gc'=>3.81, 'gg'=>4.08, 'gt'=>3.41, 'ta'=>3.21, 'tc'=>3.47, 'tg'=>3.23, 'tt'=>3.16}
  @@ROLL = {'aa'=>2.3, 'ac'=>-2, 'ag'=>0.5, 'at'=>-8.1, 'ca'=>7.4, 'cc'=>1.4, 'cg'=>6.3, 'ct'=>0.5, 'ga'=>5, 'gc'=>-0.4, 'gg'=>1.4, 'gt'=>-2, 'ta'=>8.4, 'tc'=>5, 'tg'=>7.4, 'tt'=>2.3}
  @@TILT = {'aa'=>-0.4, 'ac'=>-0.9, 'ag'=>-2.6, 'at'=>0, 'ca'=>0.6, 'cc'=>-1.1, 'cg'=>0, 'ct'=>-2.6, 'ga'=>-0.4, 'gc'=>0, 'gg'=>-1.1, 'gt'=>-0.9, 'ta'=>0, 'tc'=>-0.4, 'tg'=>0.6, 'tt'=>-0.4}
  @@TIP = {'aa'=>1.76, 'ac'=>2, 'ag'=>0.9, 'at'=>1.87, 'ca'=>-1.64, 'cc'=>0.71, 'cg'=>0.22, 'ct'=>0.9, 'ga'=>1.35, 'gc'=>2.5, 'gg'=>0.71, 'gt'=>2, 'ta'=>6.7, 'tc'=>1.35, 'tg'=>-1.64, 'tt'=>1.76}
  @@TWIST = {'aa'=>37.6, 'ac'=>35.7, 'ag'=>35.8, 'at'=>39.7, 'ca'=>32.2, 'cc'=>35.5, 'cg'=>33.9, 'ct'=>35.7, 'ga'=>38.4, 'gc'=>37.4, 'gg'=>35.5, 'gt'=>35.8, 'ta'=>34.6, 'tc'=>38.4, 'tg'=>32.2, 'tt'=>37.6}
  @@WEDGE = {'aa'=>7.2, 'ac'=>1.1, 'ag'=>8.4, 'at'=>2.6, 'ca'=>3.5, 'cc'=>2.1, 'cg'=>6.7, 'ct'=>8.4, 'ga'=>5.3, 'gc'=>5, 'gg'=>2.1, 'gt'=>1.1, 'ta'=>0.9, 'tc'=>5.3, 'tg'=>3.5, 'tt'=>7.2}
  
  
  ##
  # EXTERNAL PROGRAM METHODS
  ##
  
  # Run Orchid on the DNA sequence and return an array of the resulting values
  # NOTE: Requires Perl
  def orchid
    return Orchid.sequence(self)
  end
  
  
  ##
  # CUSTOM METHODS (not polling web @@service)
  ##
  
  # Compute a binding affinity landscape based on a PWM
  # Largely adapted from DynaPro (Morozov et al)
  def pwm_affinity(pwm)
    values = Array.new
    self.window_search(pwm.length) do |window|
      values << -Math.log(0.5*Math.exp(-pwm.score(window)) + 0.5*Math.exp(-pwm.score(window.reverse_complement)))
    end
    
    # Pad the return result with zeros to account for PWM length
    return [Array.new(pwm.length/2, 0), values, Array.new(pwm.length/2, 0)].flatten
  end
  
  # Compute nucleosome occupancy according to Desiree Tillo algorithm
  def tillo_nucleosome_model
    values = Array.new
    
    # Compute the frequencies of all 4-mers in the first window, forward and reverse-complement
    freq = Hash.new(0)
    first_window = self.subseq(1,149)
    first_window.window_search(4) { |oligo| freq[oligo] += 1 }
    first_window.reverse_complement.window_search(4) { |oligo| freq[oligo] += 1 }
      
    # Compute GC-content, propeller twist, and slide for the middle 75bp of the first window (rolling sums)
    middle = first_window.subseq(37,111)
    rolling_gc_content = middle.gc?.sum
    rolling_propeller_twist = middle.propeller_twist.sum
    rolling_slide = middle.slide.sum
    
    # Iterate over 150 bp windows
    self.window_search(150) do |window|
      # Add the last 4-mer to the frequency hash
      final_4mer = window.subseq(147,150)
      freq[final_4mer] += 1
      freq[final_4mer.reverse_complement] += 1
      
      # Add the last nucleotide of the middle 75bp to the rolling sums
      final_nuc = window.subseq(112,112)
      final_dinuc = window.subseq(111,112)
      rolling_gc_content += 1 if final_nuc == 'g' or final_nuc == 'c'
      rolling_propeller_twist += @@PROPELLER[final_dinuc]
      rolling_slide += @@SLIDE[final_dinuc]
      
      # Compute the nucleosome occupancy and store
      values << 1.67175*rolling_gc_content/75 + 0.145742*rolling_propeller_twist/75 + 1.31928*rolling_slide/75 - 0.10549*freq['aaaa'] - 0.07628*freq['aaat'] - 0.03006*freq['aagt'] - 0.05055*freq['aata'] - 0.02564*freq['aatt'] - 0.02154*freq['agaa'] - 0.03949*freq['ataa'] - 0.02354*freq['atat'] - 0.03214*freq['atta'] - 0.03314*freq['gaaa'] - 0.0334*freq['tata'] + 1.788022
        
      # Remove the first 4-mer from the frequency hash
      first_4mer = window.subseq(1,4)
      freq[first_4mer] -= 1
      freq[first_4mer.reverse_complement] -= 1
      
      # Remove the first nucleotide of the middle 75bp from the rolling sums
      first_nuc = window.subseq(37,37)
      first_dinuc = window.subseq(37,38)
      rolling_gc_content -= 1 if first_nuc == 'g' or first_nuc == 'c'
      rolling_propeller_twist -= @@PROPELLER[first_dinuc]
      rolling_slide -= @@SLIDE[first_dinuc]
    end
    
    # Pad the return result with zeros to account for window length
    return [Array.new(75, 0), values, Array.new(74, 0)].flatten
  end
  
  
  ##
  # DINUCLEOTIDE CONSTANT METHODS (basic structural properties)
  ##  

  #  We use propeller twist as a measure of helix rigidity, since the propeller twist 
  #  angles have been shown to be inversely related to rigidity of the DNA helix in crystals 
  #  (el Hassan et al. 1996). Thus, a region with high propeller twist would 
  #  mean the helix is quite rigid in this area, and similarly regions that are quite 
  #  flexible would have a low propeller twist. Propeller twist values were obtained from 
  #  cr et et al. 1997allographic data (el et al. 1996), with the exception of the TA 
  #  step, which was taken from a theoretical estimate (Gorin et al. 1995). Plots using 
  #  other sets of propeller twist dinucleotide values were very similar (data not shown). 
  #  The average propeller twist value in the entire E. coli K-12 genome is -12.63 degrees. 
  # 
  #  * Goffeau et al. The Yeast Genome Directory (1997) 387 (supplement):5-105 
  #
  #  * M.A. el Hassan and C.R. Calladine Propeller-twisting of base-pairs and the conformational 
  #    mobility of dinucleotide steps in DNA. (1996) 259:95-103 
  #
  #  * A.A. Gorin and V.B. Zhurkin and W.K. Olson B-DNA twisting correlates with base-pair 
  #    morphology. (1995) 247:34-48 
  #
  # NOTE: Returns sequence.length - 1 values (dinucleotide)
  def propeller_twist
    dinucleotide_lookup(@@PROPELLER)
  end
  
  # Slide for each dinucleotide as from Tillo averages
  def slide
    dinucleotide_lookup(@@SLIDE)
  end
  
  # Free energy change for each dinucleotide
  def free_energy
    dinucleotide_lookup(@@FREE_ENERGY)
  end
  
  # Enthalpy change for each dinucleotide
  def enthalpy
    dinucleotide_lookup(@@ENTHALPY)
  end
  
  # Entropy change for each dinucleotide
  def entropy
    dinucleotide_lookup(@@ENTROPY)
  end
  
  # Melting temperature for each dinucleotide
  def melting
    dinucleotide_lookup(@@MELTING)
  end
  
  # Clash strength for each dinucleotide
  def clash_strength
    dinucleotide_lookup(@@CLASH_STRENGTH)
  end

  # Major groove mobility for each dinucleotide
  def major_groove_mobility
    dinucleotide_lookup(@@MAJOR_GROOVE_MOBILITY)
  end
  
  # Minor groove mobility for each dinucleotide
  def minor_groove_mobility
    dinucleotide_lookup(@@MINOR_GROOVE_MOBILITY)
  end
  
  # Major groove mobility for each dinucleotide
  def major_groove_size
    dinucleotide_lookup(@@MAJOR_GROOVE_SIZE)
  end
  
  # Major groove mobility for each dinucleotide
  def minor_groove_size
    dinucleotide_lookup(@@MINOR_GROOVE_SIZE)
  end
  
  # Rise for each dinucleotide
  def rise
    dinucleotide_lookup(@@RISE)
  end
  
  # Roll for each dinucleotide
  def roll
    dinucleotide_lookup(@@ROLL)
  end
  
  # Tilt for each dinucleotide
  def tilt
    dinucleotide_lookup(@@TILT)
  end
  
  # Tip for each dinucleotide
  def tip
    dinucleotide_lookup(@@TIP)
  end
  
  # Twist for each dinucleotide
  def twist
    dinucleotide_lookup(@@TWIST)
  end
  
  # Wedge for each dinucleotide
  def wedge
    dinucleotide_lookup(@@WEDGE)
  end
  
  #  Base-stacking energies are from the dinucleotide values provided by (Ornstein et 
  #  al. 1978). The scale is in kcal/mol, and the dinucleotide values range from -3.82 
  #  kcal/mol (will melt easily) up to a maximum value of -14.59 kcal/mol (which would 
  #  require more energy to destack or melt the helix). (All 10 values are listed in the 
  #  table below.) A positive peak in base-stacking (i.e., numbers closer to 0) reflectsregions 
  #  of the helix which would de-stack or melt more readily. Conversely, minima (larger 
  #  negative numbers) in this plot would represent more stable regions of the chromosome. 
  #   
  #       Dinucleotide melting energies in kcal/mols:
  #      
  #         (GC).(GC)  -14.59
  #         (AC).(GT)  -10.51
  #         (TC).(GA)   -9.81
  #         (CG).(CG)   -9.61
  #         (GG).(CC)   -8.26
  #         (AT).(AT)   -6.57
  #         (TG).(CA)   -6.57
  #         (AG).(CT)   -6.78
  #         (AA).(TT)   -5.37
  #         (TA).(TA)   -3.82
  # 
  #    * R.L. Ornstein and R. Rein and D.L. Breen and R.D. MacElroy An optimized potential 
  #      function for the calculation of nucleic acid interaction energies. I. Base stacking 
  #      (1978) 17:2341-2360
  #
  # NOTE: Returns sequence.length - 1 values (dinucleotide)
  def stacking_energy
    dinucleotide_lookup(@@STACKING_ENERGY)
  end
  
  # Return the number of occurrences of an oligo in the given sequence
  def frequency(oligo)
    oligo.downcase!
    count = 0
    self.window_search(oligo.length) { |seq| count += 1 if seq == oligo }
    
    return count
  end
  
  
  ##
  # BASE CONTENT METHODS (complete those available already in BioRuby)
  ##
  
  #  The "G Content" of a given sequence is merely the fraction of G's in a given sequence 
  #  (Jensen et al. 1999). It can range from 0(no G's), to 1 (all G's). For a sequence 
  #  that is 50% AT content, one would expect roughly 25% G's. 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777
  def g?
    self.split(//).map { |base| base == 'g' ? 1 : 0 }
  end
  
  # The "A Content" of a given sequence is merely the fraction of A's in a given sequence 
  # (Jensen et al. 1999). It can range from 0(no A's), to 1 (all A's). For a sequence 
  # that is 50% AT content, one would expect roughly 25% A's. 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777
  def a?
    self.split(//).map { |base| base == 'a' ? 1 : 0 }
  end
  
  #  The "C Content" of a given sequence is merely the fraction of C's in a given sequence
  #  (Jensen et al. 1999). It can range from 0(no C's), to 1 (all C's). For a sequence 
  #  that is 50% AT content, one would expect roughly 25% C's. 
  # 
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes
  #    (1999) 150:773-777 
  def c?
    self.split(//).map { |base| base == 'c' ? 1 : 0 }
  end
  
  # The "T Content" of a given sequence is merely the fraction of T's in a given sequence 
  # (Jensen et al. 1999). It can range from 0(no T's), to 1 (all T's). For a sequence 
  # that is 50% AT content, one would expect roughly 25% T's. 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777 
  def t?
    self.split(//).map { |base| base == 't' ? 1 : 0 }
  end
  
  # Return the total content/percentages for the following properties on the entire sequence (average, not rolling)
  # (analogous to :gc_percent, :gc_content, :at_content, :at_percent provided in gem Bio::Sequence::NA class)
  def a_content
    return Rational(self.count('a'), self.length)
  end
  
  def a_percent
    return (100*self.a_content).to_i
  end
  
  def t_content
    return Rational(self.count('t'), self.length)
  end
  
  def t_percent
    return (100*self.t_content).to_i
  end
  
  def c_content
    return Rational(self.count('c'), self.length)
  end
  
  def c_percent
    return (100*self.c_content).to_i
  end
  
  def g_content
    return Rational(self.count('g'), self.length)
  end
  
  def g_percent
    return (100*self.g_content).to_i
  end
  
  # The percent AT is a running average of the AT content, over a given window size. 
  # Typically for a bacterial genomes of about5 Mbp, the window size is 10,000 bp. The 
  # Percent AT can range from 0 (no AT content) to 1 (100% AT). The Percent AT iscorrelated 
  # with other DNA structural features, such that AT rich regions are often more readily 
  # melted, tend to be lessflexible and more rigid, although they can also be readily 
  # compacted chromatin proteins (Pedersen et al. 2000). 
  #
  #  * A.G. Pedersen and L.J. Jensen and H.H. St\aerfeldt and S. Brunak and D.W. Ussery 
  #    A DNA structural atlas of \textitE. coli (2000) 299:907-930
  def at?
    self.split(//).map { |base| (base == 'a' or base == 't') ? 1 : 0 }
  end
  
  # Same for GC content
  def gc?
    self.split(//).map { |base| (base == 'g' or base == 'c') ? 1 : 0 }
  end
  
  # For many genomes there is a strand bias, such that one strand tends to have more 
  # G's, whilst the other strand has more C's.This GC-skew bias can be measured the number 
  # of G's minus the number of C's over a fixed length (e.g. 10,000 bp) of DNA(Jensen 
  # et al. 1999). The values can range from +1 (all G's on the examined sequence, with 
  # all C's on the other strand), to -1(the reverse case - all C's on the examined sequence, 
  # and all G's on the other strand). There is a correlation with GC-skew and the replication 
  # leading and lagging strands. 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777 
  def rolling_gc_skew
    self.split(//).map do |base| 
      case base
        when 'g' then 1
        when 'c' then -1
        else 0
      end
    end
  end
  
  # For some genomes there is also an AT strand bias, such that one strand tends to have 
  # more A's, whilst the other strand hasmore T's. This AT-skew bias is measured as the 
  # number of A's minus the number of T's over a fixed length (e.g. 10,000 bp) ofDNA 
  # (Jensen et al. 1999). The values can range from +1 (all A's on the examined sequence, 
  # with all T's on the other strand), to-1 (the reverse case - all T's on the examined 
  # sequence, and all A's on the other strand). For some genomes, there is acorrelation 
  # with AT-skew and the replication leading and lagging strands. 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777 
  def rolling_at_skew
    self.split(//).map do |base| 
      case base
        when 'a' then 1
        when 't' then -1
        else 0
      end
    end
  end

  
  ##
  # WEB SERVICE METHODS
  ##

  #  DNA curvature is calculated using the CURVATURE programme (Bolshoy et al. 1991, Shpigelman 
  #  et al. 1993). The term curved DNA here refers to DNA that is intrinsically curved 
  #  in solution and can be readily characterised by anomalous migration in acrylamide 
  #  gels. There are different models for curved DNA (Sinden et al. 1998), although the 
  #  predictions for curvature fragments largerthan a few hundred bp is essentially the 
  #  same (Haran et al. 1994). The scale is in arbitrary "Curvature units", which ranges 
  #  from 0 (e.g. no curvature) to 1.0, which is the curvature of DNA when wrapped around 
  #  the nucleosome. The scale used for this atlas ranges 3 standard deviations around 
  #  the mean. 
  #
  #   * R.R. Sinden and C.E. Pearson and V.N. Potaman and D.W. Ussery DNA: Structure and 
  #     Function (1998) 5A:1-141 
  #
  #   * E.S. Shpigelman and E.N. Trifonov and A. Bolshoy CURVATURE: Software for the Analysis 
  #     of Curved DNA. (1993) 9:435-444 
  #
  #   * T.E. Haran and J.D. Kahn and D.M. Crothers Sequences elements responsible for 
  #     DNA curvature (1994) 225:729-738 
  #
  #   * A. Bolshoy and P. McNamara and R.E. Harrington and E.N. Trifonov Curved DNA Without 
  #     A-A - Experimental Estimation of All 16 DNA Wedge Angles (1991) 88:2312-2316 
  #
  #  NOTE: Returns sequence.length - 20 values, pad with zeros to account
  def curvature
    [Array.new(10, 0), compute('Intrinsic Curvature'), Array.new(10, 0)].flatten
  end
  
  # A trinucleotide model based on the preferential location 
  # of sequences within nucleosomal core sequences (Satchwell et al. 1986). We use the 
  # magnitude (e.g.absolute values) of the trinucleotide numbers as a measure of DNA 
  # flexibility (Baldi et al. 1996). The trinucleotide values range from essentially 
  # zero (0.003, presumably more flexible), to 0.28 (considered rigid). Since very few 
  # of the trinucleotide have values close to zero (e.g. little preference for nucleosome 
  # positioning), this measureis considered most sensitive towards the low ("flexibity") 
  # 
  #
  #  * S.C. Satchwell and H.R. Drew and A.A. Travers Sequence periodicities in chicken 
  #    nucleosome core DNA (1986) 191:659-675 
  #
  #  * P. Baldi and S. Brunak and Y. Chauvin and A. Krogh Naturally occurring nucleosome 
  #    positioning signals in human exons and introns. (1996) 263:503-510 
  #
  # NOTE: Returns sequence.length - 2 values, pad with zeros to account
  def nucleosome_preference
    [0, compute('Position Preference'), 0].flatten
  end
  
  # "Protein Induced Deformability" dinucleotide values are from protein induced deformation 
  #  of DNA helices as determined by examination of more than a hundred cr et et al. 1997al 
  #  structures of DNA/protein complexes (Olson et al. 1998). The dinucleotide values 
  #  range from 2.1 (the least deformable dinucleotide), to 12.1 (i.e., the dinucleotide 
  #  step (CpG), which is often deformed by proteins). Thus, on this scale, a larger value 
  #  reflects a more deformable sequence whilst a smaller value indicates a region where 
  #  the DNA helix is less likely to be changed dramatically by proteins. The average 
  #  protein deformability value in the entire E. coli K-12 genome is 5.12. 
  #
  #  * Goffeau et al. The Yeast Genome Directory (1997) 387 (supplement):5-105 
  #
  #  * W.K. Olson and A.A. Gorin and X.J. Lu and L.M. Hock and V.B. Zhurkin DNA sequence-dependent 
  #    deformability deduced from protein-DNA crystal complexes. (1998) 95:11163-11168 
  #
  # NOTE: Returns sequence.length - 1 values
  def protein_deformability
    compute('Protein Deformability')
  end
  
  #  DNase I values are based on experimentally determined trinucleotide values (Brukner 
  #  et al. 1995, Brukner et al. 1995). These values are reflective of the anisotropic 
  #  flexibility or "bendability" of a particular DNAsequence. The trinucleotide values 
  #  range from -0.280 (rigid) to +0.194 (very "bendable" towards the major groove). Smoothing 
  #  over a large regions, (which is necessary for viewing entire genomes) tends to smooth 
  #  out differences in bendability. The average DNase I ("bendability") value in the 
  # 
  #  * I. Brukner and R. Sanchez and D. Suck and S. Pongor Sequence-dependent bending 
  #    propensity of DNA as revealed by DNase I: parameters for trinucleotides. (1995) 14:1812-1818 
  #                                
  #
  #  * I. Brukner and R. Sanchez and D. Suck and S. Pongor Trinucleotide models for DNA 
  #    bending propensity: comparison of models based on DNaseI digestion and nucleosome 
  #    packaging data. (1995) 13:309-317
  #
  # NOTE: Returns sequence.length - 2 values, pad with zeros to account
  def dnase_sensitivity
    [0, compute('DNase I Sensitivity'), 0].flatten
  end
  
  #  For a given sequence, any palindrome of 6 nt (e.g., AAATTT) is given a value of 1, 
  #  while all bases not included in palindromic hexamers are given a value of 0 (van et 
  #  al. 2003). 
  #
  #  * van Noort V, Worning P, Ussery DW, Rosche WA, Sinden RR Strand misalignments lead 
  #    to quasipalindrome correction (2003) 19:365-9 
  def palindromic_hexamers
    compute('Palindromic hexamers')
  end
  
  # Global Direct repeats are found by taking the first 100 bp of sequence, and
  # looking for the best match within the whole segment, on the same strand, in the
  # same direction [5' to 3'] (Skovgaard et al. 2002). Values are binned into 10
  # values, and represent the lower end of the best match, and range from 0 (10% or
  # less match) to 9 (at least 90 out of the 100 nucleotides match perfectly).
  def global_direct_repeats
    compute('Global Direct Repeats')
  end
  
  # Global Direct repeats are found by taking the first 100 bp of sequence, and
  # looking for the best match within the whole segment, on the opposite strand, in
  # the same direction  [5' to 3'] (Skovgaard et al. 2002). Values are binned into
  # 10 values, and represent the lower end of the best match and range from 0 (10%
  # or less match) to 9 (at least 90 out of the 100 nucleotides match perfectly).
  #
  # * M. Skovgaard and L.J. Jensen and C. Friis and H.H. Staerfeldt,and P. Worning
  # and S. Brunak The Atlas Visualization of Genomewide Information (2002) 33:49-63
  def global_inverted_repeats
    compute('Global Inverted Repeats')
  end
  
  # Local Direct repeats are found by taking a 100 bp sequence window, and looking for 
  # the best match of a 30 bp piece withinthat window, on the same strand, in the same 
  # direction (Jensen et al. 1999). Values can range from 0 (no match at all) to 1(one 
  # or more perfect match within the window). 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777 
  def direct_repeats
    compute('Direct Repeats')
  end
  
  # Local Everted repeats are found by taking a 100 bp sequence window, and looking for 
  # the best match of a 30 bp piece withinthat window, on the opposite strand, in the 
  # same direction (Jensen et al. 1999). Values can range from 0 (no match at all) to 
  # 1(one or more perfect match within the window). 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777
  def everted_repeats
    compute('Everted Repeats')
  end
  
  # Local Inverted repeats are found by taking a 100 bp sequence window, and looking 
  # for the best match of a 30 bp piece withinthat window, on the opposite strand, in 
  # the opposite direction (Jensen et al. 1999). Values can range from 0 (no match at 
  # all)to 1 (one or more perfect match within the window). 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777
  def local_inverted_repeats
    compute('Local Inverted Repeats')
  end
  
  # Local Mirror repeats are found by taking a 100 bp sequence window, and looking for 
  # the best match of a 30 bp piece withinthat window, on the same strand, in the opposite 
  # direction (Jensen et al. 1999). Values can range from 0 (no match at all) to 1(one 
  # or more perfect match within the window). 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777
  def mirror_repeats
    compute('Mirror Repeats')
  end
  
  # "Quasi-palindromes" are short inverted repeats, which are found by taking a 30 bp 
  # piece of sequence, and looking for matcheswith at least 6 out of 7 nt matching, on 
  # the opposite strand, in the opposite direction (van et al. 2003). Values canrange 
  # from 0 (no match at all) to 1 (one or more perfect match within the window). 
  #
  #  * van Noort V, Worning P, Ussery DW, Rosche WA, Sinden RR Strand misalignments lead 
  #    to quasipalindrome correction (2003) 19:365-9 
  def quasi_palindromes
    compute('Quasi-palindromes')
  end
  
  # "Perfect-palindromes" are short inverted repeats, which are found by taking a 30 
  # bp piece of sequence, and looking forperfect matches of 7 nt or longer, on the opposite 
  # strand, in the opposite direction (van et al. 2003). Values can rangefrom 0 (no match 
  # at all) to 1 (one or more perfect match within the window). 
  #
  #  * van Noort V, Worning P, Ussery DW, Rosche WA, Sinden RR Strand misalignments lead 
  #    to quasipalindrome correction (2003) 19:365-9
  def perfect_palindromes
    compute('Perfect-palindromes')
  end
  
  # A "simple repeat" is a region which contains a simple oligonucleotide repeat, like 
  # microsattelites. Simple repeats are foundby looking for tandem repeats of length 
  # R within a 2R-bp window. By using the values 12, 14, 15, 16, and 18 for R, allsimple 
  # repeats of lengths 1 through 9 are calculated, of length of at least 24 bp (Jensen 
  # et al. 1999). Values can range from 0(no match at all) to 1 (one or more perfect 
  # match within the window). 
  #
  #  * L. J. Jensen and C. Friis and D.W. Ussery Three views of complete chromosomes 
  #    (1999) 150:773-777 
  def simple_repeats
    compute('Simple Repeats')
  end
  
  # Current undocumented/unimplemented web service properties are:
  #      AAAA
  #      CCCC
  #      TTTT
  #      GGGG
  #      T4 or C4 vs. A4 or G4
  #      (Y)10 vs. (R)10
  #      (CR)5 vs. (YG)5
  #      (CA)3
  #      (CG)3
  #      (TA)3
  #      (TG)3
  #      (YR)5
  # Current properties available via webservice but unused/implemented locally are:
  #      A/C/G/T Content
  #      Percent AT
  #      AT/GC Skew
  #      Propeller Twist
  #      Stacking Energy
  
  
  ##
  # Helper methods to submit a DNA Property request to the queue
  ##
  private

  # Submit a request for a specific DNA property on a sequence
  def compute(method)
    fetch_result(submit_request(method).queueentry.jobid).output.values.split(',').map { |elem| elem.to_f }
  end 
  
  def submit_request(method)
    @@service.DNApropertyRun(:parameters => {:method => method, :sequence => {:id => 'NoID', :seq => self}})
  end
  
  # Fetches the result of a job placed in the queue
  def fetch_result(jobid)
    sleep 5 until poll_queue(jobid) == 'FINISHED'
    @@service.DNApropertyFetchResult(:job => {:jobid => jobid})
  end
  
  # Polls the server queue to see if a job is finished
  # Possible values are QUEUED, ACTIVE, FINISHED, WAITING, REJECTED, UNKNOWN JOBID or QUEUE DOWN
  def poll_queue(jobid)
    @@service.pollQueue(:job => {:jobid => jobid}).queueentry.status
  end
  
  # Map dinucleotides in a hash and return the result array
  def dinucleotide_lookup(lookup_table)
    values = Array.new
    self.window_search(2) do |dinucleotide|
      values << lookup_table[dinucleotide]
    end
    
    return values
  end
end