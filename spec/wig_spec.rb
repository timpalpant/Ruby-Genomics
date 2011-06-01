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

describe WigFile do

  before do
    @test = WigFile.new(TEST_WIG)
  end
  
  it "should find all chromosomes" do
    @test.chromosomes.length.should == 3
    @test.include?('chrXI').should be_true
    @test.include?('chr1').should be_true
    @test.include?('2micron').should be_true
  end
  
  it "should correctly index all chromosome starts" do
    @test.chr_start('chrXI').should == 2
    @test.chr_start('chr1').should == 29
    @test.chr_start('2micron').should == 45
  end
  
  it "should correctly index all chromosome stops" do
    @test.chr_stop('chrXI').should == 28
    @test.chr_stop('chr1').should == 44
    @test.chr_stop('2micron').should == 50
  end
  
  it "should correctly index all chromosome data lengths" do
    @test.chr_length('chrXI').should == 26
    @test.chr_length('chr1').should == 15
    @test.chr_length('2micron').should == 5
  end
  
  it "should load whole fixedStep chromosomes" do
    data = @test['chrXI']
    data.length.should == 26
    
    data = @test.chr('chr1')
    data.length.should == 15
  end
  
  it "should load whole variableStep chromosomes" do
    data = @test['2micron']
    data.length.should == 12
  end
  
  it "should iterate over all chromosomes" do
    count = 0
    @test.each { |chr, values| count += 1 }
    count.should == 3
  end
  
  it "should query randomly within fixedStep chromosomes" do
    @test.query('chr1', 5, 8).should == [5, 6, 7, 8]
    lambda { @test.query('chrXI', 25, 35) }.should raise_error
  end
  
  it "should raise an error when attempting to query a variableStep chromosome" do
    lambda { @test.query('2micron', 101, 110) }.should raise_error
  end

end