#
#  read_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'read'

describe Read do
  before do
    @test = Read.new('chrI', 1, 36)
  end
  
  it "should correctly output to SAM format"
  
  it "should correctly output to Bowtie hits format"
end