#
#  assembly_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'assembly'

describe Assembly do
  before do
    @yeast = Assembly.load('sacCer2')
  end
  
  it "should load all chromosomes" do
    @yeast.chromosomes.length.should == 18
  end
  
  it "should have the correct length for each chromosome" do
    @yeast['chrIV'].should ==	1531919
    @yeast['chrXV'].should ==	1091289
    @yeast['chrVII'].should ==	1090947
    @yeast['chrXII'].should ==	1078175
    @yeast['chrXVI'].should ==	948062
    @yeast['chrXIII'].should ==	924429
    @yeast['chrII'].should ==	813178
    @yeast['chrXIV'].should ==	784333
    @yeast['chrX'].should ==	745742
    @yeast['chrXI'].should ==	666454
    @yeast['chrV'].should ==	576869
    @yeast['chrVIII'].should ==	562643
    @yeast['chrIX'].should ==	439885
    @yeast['chrIII'].should ==	316617
    @yeast['chrVI'].should ==	270148
    @yeast['chrI'].should ==	230208
    @yeast['chrM'].should ==	85779
    @yeast['2micron'].should ==	6318
  end
  
  it "should return its len file" do
    @yeast.len_file.should == File.expand_path(File.dirname(__FILE__) + "/../resources/assemblies/sacCer2.len")
  end
end