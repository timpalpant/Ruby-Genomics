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

end

module BGZip

end