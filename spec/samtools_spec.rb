#
#  samtools_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 5/31/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'samtools'

TEST_BAM = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bam')

describe SAMTools do
	it "should have the correct number alignments for each chromosome"
end