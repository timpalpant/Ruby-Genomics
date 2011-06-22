require 'enumerator'

##
# A line-oriented data file
# base class for Bed, BedGraph, GFF, SAM, etc.
# Idea is that it should behave somewhat like a File object,
# but that relevant parsed Entry objects (BedEntry, SAMEntry, etc.) 
# are returned instead of lines
##
class EntryFile
  include Enumerable

  def initialize(filename)
    @data_file = File.expand_path(filename)
  end
  
  # Perform any additional cleanup operations (deleting indexes, etc.)
  def close
  end
  
  # Open the EntryFile (optionally with a block)
  def self.open(filename, &block)
    entry_file = self.new(filename)
    
    if block
      yield entry_file
      entry_file.close
    else
      return entry_file
    end
  end
  
  # Iterate over each of the entries in an EntryFile
  def self.foreach(filename, chr = nil, start = nil, stop = nil)
    entry_file = self.new(filename)
    entry_file.each(chr, start, stop) { |entry| yield entry }
    entry_file.close
  end
  
  # Make entry files enumerable
  def each(chr = nil, start = nil, stop = nil)
    query_lines(chr, start, stop) do |line|
      # Skip comment and track lines
      next if line.start_with?('#') or line.start_with?('@') or line.start_with?('track') or line.chomp.empty?
      yield parse(line)
    end
  end
  
  def to_bed(output)
    to_disk(output) { |entry| entry.to_bed }
  end
  
  def to_bedgraph(output)
    to_disk(output) { |entry| entry.to_bedgraph }
  end
  
  def to_gff(output)
    to_disk(output) { |entry| entry.to_gff }
  end
  
  private
  
  # Should be overridden in subclasses to parse an line into an object
  def parse(line)
    raise EntryFileError, "Do not know how to parse the entries in #{File.basename(@data_file)}!"
  end
  
  # Should be overridden in subclasses to query the entry file for lines
  def query_lines(chr = nil, start = nil, stop = nil)
    raise EntryFileError, "Do not know how to query for lines in #{File.basename(@data_file)}!"
  end
  
  def to_disk(output)
    File.open(File.expand_path(output), 'w') do |f|
      self.each do |entry|
        f.puts yield(entry)
      end
    end
  end
end

class TextEntryFile < EntryFile
  private
  
  # Get all lines in the file matching chr:start-stop
  def query_lines(chr = nil, start = nil, stop = nil)
    raise "Querying TextEntryFiles for specific intervals (start-stop) is not yet implemented" if start or stop
    
    if chr.nil?
      File.foreach(@data_file) { |line| yield line }
    else
      IO.popen("grep -w #{chr} #{File.expand_path(filename)}") do |output|
        output.each { |line| yield line }
      end
    end
  end
end

class BinaryEntryFile < EntryFile
  private
  
  # Query the binary file and return the resulting text-entry lines
  def query_lines(chr, start, stop)
    IO.popen(query_command(chr, start, stop)) do |output|
      output.each { |line| yield line }
    end
  end
  
  # Should be overridden in subclasses to construct the query command
  def query_command(chr = nil, start = nil, stop = nil)
    raise "Do not know how to query binary file #{File.basename(@data_file)}"
  end
end

class EntryFileError < StandardError
end