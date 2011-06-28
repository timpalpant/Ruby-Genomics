#
#  wig_spec.rb
#  BioRuby
#
#  Created by Timothy Palpant on 4/8/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'wig'

TEST_WIG = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.wig')
TEST_BIGWIG = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bw')

# test.wig and test.bw have the same contents
# so many of the tests can be shared
shared_examples "wig file" do
  context "when indexing" do
    it "should find all chromosomes" do
      @test.chromosomes.length.should == 3
      @test.include?('chrXI').should be_true
      @test.include?('chrI').should be_true
      @test.include?('2micron').should be_true
    end
    
    it "should index all chromosome starts" do
      @test.chr_start('chrXI').should == 20
      @test.chr_start('chrI').should == 1
      @test.chr_start('2micron').should == 100
    end
    
    it "should index all chromosome stops" do
      @test.chr_stop('chrXI').should == 148
      @test.chr_stop('chrI').should == 15
      @test.chr_stop('2micron').should == 111
    end
  end
  
  context "when iterating" do
    it "should iterate over all contigs" do
      count = 0
      @test.each { |contig| count += 1 }
      count.should == 3
    end
    
    it "should iterate over all chunks" do
      count = 0
      @test.each { |chunk| count += 1 }
      count.should == 3
    end
  end
  
  context "when querying" do
    it "should query randomly within fixedStep chromosomes" do
      result = @test.query('chrI', 5, 8)
      result.length.should == 4
      (5..8).each { |bp| result[bp].should == bp }
      
      result = @test.query('chrXI', 25, 35)
      result.length.should == 11
      (25..28).each { |bp| result[bp].should == 3 }
      result[29].should be_nil
      (30..33).each { |bp| result[bp].should == 4 }
      result[34].should be_nil
      result[35].should == 9
    end
    
    it "should query randomly within variableStep chromosomes" do
      result = @test.query('2micron', 101, 110)
      result.length.should == 10
      result[101].should == 6
      result[102].should be_nil
      result[103].should be_nil
      result[104].should be_nil
      result[105].should == 10
      result[106].should be_nil
      result[107].should be_nil
      result[108].should be_nil
      result[109].should be_nil
      result[110].should == 1
    end
  end
end

describe BigWigFile do
  before do
    @test = BigWigFile.new(TEST_BIGWIG)
  end
  
  after do
    @test.close
  end
  
  include_examples "wig file"
  
  it "should allow opening with a block" do
    BigWigFile.open(TEST_BIGWIG) { |bw| }
  end
  
  it "should allow querying for the average value of a window" do
    @test.query_average('2micron', 100, 101).should == 5.5
    @test.query_average('chrI', 1, 3).should == 2
  end
  
  it "should output to BedGraph"
  it "should output to Wig"
end

describe WigFile do
  before do
    @test = WigFile.new(TEST_WIG)
  end
  
  after do
    @test.close
  end
  
  include_examples "wig file"
  
  it "should allow opening with a block" do
    WigFile.open(TEST_WIG) { |wig| }
  end
  
  it "should output to BedGraph"
  it "should output to BigWig"
end