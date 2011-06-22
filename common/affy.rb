require 'entry_file'
require 'spot_file'
require 'spot'

# Affy file entries are just spots
class AffyEntry < Spot
end

##
# Load Affymetrix array files
##
class AffyFile < TextEntryFile
  extend SpotFile

  # Override each because Affy files are not strictly line-based
  # We have to look at multiple lines to parse individual entries
  def each(query_chr = nil, query_start = nil, query_stop = nil)
		prev_chr, prev_start, prev_value = nil, nil, nil
		File.foreach(@data_file) do |line|
			# Skip comment and empty lines
			next if line.start_with?('#') or line.chomp.empty?
      
			entry = line.chomp.split("\t")
      raise AffyFileError, "Not a valid Affy file entry!" if entry.length < 3
      
			chr = entry[0]
			start = entry[1].to_i
			value = entry[2].to_f
    
			# Extend values (spots) up to the next spot (i.e. tile)
			if prev_chr and prev_start and prev_value and prev_chr == chr and prev_start < start
				spot = AffyEntry.new
        spot.chr = prev_chr
				spot.start = prev_start
				spot.stop = start-1
				spot.value = prev_value
				
        if query_chr.nil?
          yield spot
        else
          yield spot if spot.chr == query_chr
        end
			end
		
			prev_chr = chr
			prev_start = start
			prev_value = value
		end
  end
end

class AffyFileError < StandardError
end