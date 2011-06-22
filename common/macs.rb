require 'entry_file'
require 'spot_file'
require 'spot'

##
# An entry in a MACS output file
##
class MACSEntry < Spot
	attr_accessor :summit, :tags, :pvalue, :fdr
	
	def self.parse(line)
    entry = line.chomp.split("\t")

    spot = self.new
		spot.chr = entry[0]
    spot.start = entry[1].to_i
    spot.stop = entry[2].to_i
    spot.summit = spot.start + entry[4].to_i
    spot.tags = entry[5].to_i
    spot.pvalue = 10**(entry[6].to_f / -10)
    spot.value = entry[7].to_f
    spot.fdr = entry[8].to_f
      
    return spot
	end
  
  def fold_enrichment
    @value
  end
end

class MACSFile < TextEntryFile
  extend SpotFile

  private
  
  def parse(line)
    MACSEntry.parse(line)
  end
end
