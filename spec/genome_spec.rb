require 'spec_helper'
require 'genome'

describe Genome do
  TEST_2BIT = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.2bit')
  
  before do
    @test = Genome.new(TEST_2BIT)
  end
  
  it "should index all chromosomes" do
    @test.length.should == 2
  end
  
  it "should index all chromosome lengths" do
    @test['chrI'].should == 999
    @test['chrXI'].should == 4999
  end
  
  it "should include the chromosome chrII" do
    @test.include?('chrXI').should be_true
  end
  
  it "should not include the chromosome chr8" do
    @test.include?('chr8').should be_false
  end
  
  it "should include the interval chrI:1-100" do
    @test.include?('chrI', 1, 100).should be_true
  end
  
  it "should not include the interval chrI:1-20,000" do
    @test.include?('chrI', 1, 20_000).should be_false
  end
  
  it "should be able to query for sequences" do
    seq = @test.sequence('chrI', 10, 20)
    seq.should == Bio::Sequence::NA.new('CCCACACACCC')
  end
  
  it "should raise an error if the sequences are not available" do
    lambda { @test.sequence('chr8', 1, 20) }.should raise_error
  end
  
  it "should count the total number of base pairs" do
    @test.num_bases.should == 5998
  end
end
