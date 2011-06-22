##
# A file with Read entries
# Extended in EntryFile subclasses
# Don't allow loading all reads, but allow iterating specifically over read pairs
##
module ReadFile
  # Return each read, but only one (forward) entry for paired-end reads
  def foreach_read(filename, chr = nil, start = nil, stop = nil)
    self.foreach(filename, chr, start, stop) do |entry|
      yield entry unless entry.paired? and entry.crick?
    end
  end
end