#
#  affy_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 5/31/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'affy'

describe AffyFile do
  AFFY_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.affy')
  
  before do
    @test = AffyFile.new(AFFY_FILE)
  end
  
  after do
    @test.close
  end
  
  it "should allow opening with a block" do
    AffyFile.open(AFFY_FILE) { |bed| }
  end
  
  it "should allow opening without a block" do
    affy = AffyFile.open(AFFY_FILE)
    affy.close
  end
  
  context "when iterating over all entries" do
    it "should have 10 entries" do
      @test.count.should == 10
    end
    
    it "should return the number of skipped entries"
    
    it "should have 4 chromosomes" do
      @test.chromosomes.length.should == 4
    end
    
    it "should iterate over all the entries" do
      count = 0
      @test.each { |entry| count += 1 }
      count.should == 10
    end
    
    it "should allow static iteration over all the entries" do
      count = 0
      AffyFile.foreach(AFFY_FILE) { |entry| count += 1 }
      count.should == 10
    end
  end
  
  context "when querying for a specific chromosome" do
    it "should have 2 entries on chrI" do
      @test.count('chrI').should == 2
      
      count = 0
      @test.each('chrI') { |entry| count += 1 }
      count.should == 2
      
      @test.chr('chrI').length.should == 2
      @test['chrI'].length.should == 2
      
      count = 0
      @test.chr('chrI') { |entry| count += 1 }
      count.should == 2
    end
    
    it "should have 0 entries on chr8" do
      @test.count('chr8').should == 0
      
      count = 0
      @test.each('chr8') { |entry| count += 1 }
      count.should == 0
      
      @test.chr('chr8').length.should == 0
      @test['chr8'].length.should == 0
      
      count = 0
      @test.chr('chr8') { |entry| count += 1 }
      count.should == 0
    end
  end
  
  context "when querying for a specific interval" do
    it "should count the correct number of entries" do
      @test.count('chrIV', 4, 9).should == 2
    end
    
    it "should return the correct number of entries" do      
      count = 0
      @test.each('chrIV', 4, 9) { |entry| count += 1 }
      count.should == 2
    end
    
    it "should return the correct entries" do
      @test.each('chrIV', 4, 9) do |entry|
        entry.chr.should == 'chrIV'
        entry.high.should >= 4
        entry.low.should <= 9
      end
    end
  end
  
  context "when querying for a specific spot id" do
    it "should raise an error" do
      lambda { @test.id('MyAwesomeSpot') }.should raise_error
    end
  end
  
  context "when computing statistics" do
    it "should have total = 49.8" do
      @test.total.should == 49.800000000000004
    end
    
    it "should have mean = 4.98" do
      @test.mean.should == 4.98
    end
    
    it "should have standard deviation = 3.9773860763069" do
      @test.stdev.should == 3.9773860763068996
    end
  end
  
  it "should allow querying for values" do
    result = @test.query('chrIII', 10, 18)
    result.start.should == 15
    result.stop.should == 18
    result.length.should == 4
    result.coverage.should == 4
    (10..14).each { |bp| result.get(bp).should be_nil }
    (15..18).each { |bp| result[bp].should == 2.0 }
  end
end
