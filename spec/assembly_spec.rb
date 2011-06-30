#
#  assembly_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'assembly'

shared_examples "sacCer2 assembly" do
  it "should load all chromosomes" do
    @test.chromosomes.length.should == 18
  end
  
  it "should have the correct length for each chromosome" do
    @test['chrIV'].should ==	1531919
    @test['chrXV'].should ==	1091289
    @test['chrVII'].should ==	1090947
    @test['chrXII'].should ==	1078175
    @test['chrXVI'].should ==	948062
    @test['chrXIII'].should ==	924429
    @test['chrII'].should ==	813178
    @test['chrXIV'].should ==	784333
    @test['chrX'].should ==	745742
    @test['chrXI'].should ==	666454
    @test['chrV'].should ==	576869
    @test['chrVIII'].should ==	562643
    @test['chrIX'].should ==	439885
    @test['chrIII'].should ==	316617
    @test['chrVI'].should ==	270148
    @test['chrI'].should ==	230208
    @test['chrM'].should ==	85779
    @test['2micron'].should ==	6318
  end
end

describe Assembly do
  TEST_LEN_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.len')

  context "using a built-in assembly file" do
    before do
      @test = Assembly.load('sacCer2')
    end
    
    it "should return its len file" do
      @test.len_file.should == File.expand_path(File.dirname(__FILE__) + "/../resources/assemblies/sacCer2.len")
    end
    
    include_examples "sacCer2 assembly"
  end
  
  context "using a specified assembly file" do
    before do
      @test = Assembly.load(TEST_LEN_FILE)
    end
    
    it "should return its len file" do
      @test.len_file.should == TEST_LEN_FILE
    end
    
    include_examples "sacCer2 assembly"
  end
end
