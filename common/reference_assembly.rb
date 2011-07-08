#
#  reference_assembly.rb
#  BioRuby
#  Helper to check if assemblies are built-in or specified by file
#
#  Created by Timothy Palpant on 4/8/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'bio/genomics/assembly'

module ReferenceAssembly
  BUILT_IN_ASSEMBLIES_DIR = File.dirname(__FILE__) + '/../resources/assemblies'

  # Look for the desired assembly in the built-in resources directory
  # or in a specified file
  def self.load(a)
    builtin = BUILT_IN_ASSEMBLIES_DIR + "/#{a}.len"
    if File.exist?(builtin)
      return Bio::Genomics::Assembly.load(builtin)
    elsif File.exist?(File.expand_path(a))
      begin
        Bio::Genomics::Assembly.load(File.expand_path(a))
      rescue
        raise ReferenceAssemblyError, "Error loading assembly from file #{File.basename(a)}!"
      end
    else
      raise ReferenceAssemblyError, "Specified reference assembly (#{a}) is neither built-in nor a file! If you wish to add your assembly, place the len file in resources/assemblies/*"
    end
  end
end

class ReferenceAssemblyError < StandardError
end
