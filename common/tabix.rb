#
#  tabix.rb
#  ruby-genomics
#  Wrapper for the tabix and bgzip executables
#
#  Created by Timothy Palpant on 5/30/11.
#  Copyright 2011 UNC. All rights reserved.
#

# For documentation, see: http://samtools.sourceforge.net/tabix.shtml
# Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file in.tab.bgz.tbi 
# when region is absent from the command-line. The input data file must be position sorted and compressed 
# by bgzip which has a gzip(1) like interface. After indexing, tabix is able to quickly retrieve data 
# lines overlapping regions specified in the format "chr:beginPos-endPos". Fast data retrieval also works 
# over network if URI is given as a file name and in this case the index file will be downloaded if it is 
# not present locally.

module Tabix
  def self.index(filename, chr_col = 1, start_col = 2, end_col = 3)
    %x[ tabix -s #{chr_col} -b #{start_col} -e #{end_col} #{filename} ]
  end
  
  def self.query(filename, chr, start, stop, &block)
    if block
      IO.popen("tabix #{filename} #{chr}:#{start}-#{stop}") do |output|
        output.each { |line| yield line }
      end
    else
      %x[ tabix #{filename} #{chr}:#{start}-#{stop} ].split("\n")
    end
  end
end

module BGZip
  def self.compress(input_file, output_file)
    %x[ bgzip -c #{input_file} > #{output_file} ]
  end
  
  def self.decompress(input_file, output_file)
    %x[ bgzip -c -d #{input_file} > #{output_file} ]
  end
end