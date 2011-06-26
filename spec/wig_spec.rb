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

shared_examples "wig files" do
  
end

describe BigWigFile do
  
end

describe WigFile do

  before do
    @test = WigFile.new(TEST_WIG)
  end
  
  after do
    @test.close
  end
  
  it "should find all chromosomes" do
    @test.chromosomes.length.should == 3
    @test.include?('chrXI').should be_true
    @test.include?('chr1').should be_true
    @test.include?('2micron').should be_true
  end
  
  it "should correctly index all chromosome starts" do
    @test.chr_start('chrXI').should == 20
    @test.chr_start('chr1').should == 1
    @test.chr_start('2micron').should == 100
  end
  
  it "should correctly index all chromosome stops" do
    @test.chr_stop('chrXI').should == 148
    @test.chr_stop('chr1').should == 15
    @test.chr_stop('2micron').should == 111
  end
  
  it "should query randomly within fixedStep chromosomes" do
    result = @test.query('chr1', 5, 8)
    (5..8).each { |bp| result[bp].should == bp }
    
    result = @test.query('chrXI', 25, 35)
    (25..28).each { |bp| result[bp].should == 3 }
    result[29].should be_nil
    (30..33).each { |bp| result[bp].should == 4 }
    result[34].should be_nil
    result[35].should == 9
  end
  
  it "should query randomly within variableStep chromosomes" do
    result = @test.query('2micron', 101, 110)
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

describe ContigInfo do
  before do
    @test = ContigInfo.new()
  end
end