require 'spot_array'

##
# A file with Spot entries (Genomic Interval + Value), i.e. microarray
# Extended in EntryFile subclasses
##
module SpotFile
  def load(filename, chr = nil, start = nil, stop = nil)
    spot_array = SpotArray.new
  
    self.foreach(filename, chr, start, stop) do |entry|
      spot_array[entry.chr] ||= Array.new
      spot_array[entry.chr] << entry
    end
    
    return spot_array
  end
end