require 'enumerator'
require 'unix_file_utils'
require 'tabix'

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
    skipped = 0
    
    query_lines(chr, start, stop) do |line|
      # Skip comment and track lines
      next if line.start_with?('#') or line.start_with?('@') or line.start_with?('track') or line.chomp.empty?
      
      begin
        yield parse(line)
      rescue EntryFileError
        skipped += 1
      end
    end
    
    puts "Skipped #{skipped} invalid entries" if ENV['DEBUG']
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
  def initialize(filename, chr_col, start_col, end_col)
    super(filename)
    
    @chr_col = chr_col
    @start_col = start_col
    @end_col = end_col
    
    @sorted_file = @data_file + '.sorted'
    @bgzipped_file = @data_file + '.bgz'
    @index_file = @bgzipped_file + '.tbi'
  end
  
  # Delete the Tabix index file and the BGZipped version, if it exists
  def close
    File.delete(@index_file) if File.exist?(@index_file)
    File.delete(@bgzipped_file) if File.exist?(@bgzipped_file)
  end
  
  private
  
  # Get all lines in the file matching chr:start-stop
  def query_lines(chr = nil, start = nil, stop = nil)
    raise EntryFileError, "Tabix only supports queries with start AND stop" if start and stop.nil?
    
    # If we're getting all entries, just use File#foreach
    if chr.nil?
      File.foreach(@data_file) { |line| yield line }
    # If we're getting a specific chromosome, use grep to filter the entries
    elsif start.nil? or stop.nil?
      IO.popen("grep -w #{chr} #{@data_file}") do |output|
        output.each { |line| yield line }
      end
    # If we're querying for a specific region, use Tabix to index the file
    else
      index() if not File.exist?(@index_file)
      IO.popen("tabix #{@bgzipped_file} #{chr}:#{start}-#{stop}") do |output|
        output.each { |line| yield line }
      end
    end
  end
  
  # Index all TextEntryFiles with Tabix
  def index
    # File must be sorted
    File.sort(@data_file, @sorted_file, "-k#{chr_col},#{chr_col} -k#{start_col},#{start_col}n")
    
    # and BGZipped
    BGZip.compress(@sorted_file, @bgzipped_file)
    # Delete the temporary sorted text file
    File.delete(@sorted_file)
    
    # Now Tabix can index it
    Tabix.index(@bgzipped_file, @chr_col, @start_col, @end_col)
  end
end

class BinaryEntryFile < EntryFile
  def initialize(filename, index_file)
    super(filename)
    @index_file = File.expand_path(index_file)
  end
  
  # Cleanup the index, if it exists
  def close
    File.delete(@index_file) if File.exist?(@index_file)
  end
  
  private
  
  # Query the binary file and return the resulting text-entry lines
  def query_lines(chr, start, stop)
    index() if not File.exist?(@index_file)
    
    IO.popen(query_command(chr, start, stop)) do |output|
      output.each { |line| yield line }
    end
  end
  
  # Should be overridden in subclasses to construct the query command
  def query_command(chr = nil, start = nil, stop = nil)
    raise "Do not know how to query binary file #{File.basename(@data_file)}"
  end
  
  # Should be overridden in subclasses
  def index
    raise EntryFileError, "Do not know how to index binary file #{File.basename(@data_file)}"
  end
end

class EntryFileError < StandardError
end