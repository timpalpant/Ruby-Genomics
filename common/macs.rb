require 'entry_file'
require 'spot_file'
require 'spot'

##
# An entry in a MACS output file
##
class MACSEntry < Spot
  attr_accessor :summit, :tags, :pvalue, :fdr
  
  def self.parse(line)
    begin
      entry = line.chomp.split("\t")
      raise MACSEntryError, "Invalid MACS Entry: MACS peak file must have 9 columns" if entry.length != 9
      raise MACSEntryError, "MACS header detected!" if entry[0] == 'chr' and entry[1] == 'start' and entry[2] == 'end'

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
    rescue
      raise MACSEntryError, "Invalid MACS entry!"
    end
  end
  
  def fold_enrichment
    @value
  end
  
  # Can't query by id in MACS peak files since they have no id
  def id(query_id)
    raise MACSEntryError, "Cannot query for a MACS peak by id!"
  end
end

class MACSFile < TextEntryFile  
  CHR_COL = 1
  START_COL = 2
  END_COL = 3
  
  def initialize(filename)
    super(filename, CHR_COL, START_COL, END_COL)
  end

  private
  
  def parse(line)
    MACSEntry.parse(line)
  end
end

class MACSEntryError < EntryFileError
end
