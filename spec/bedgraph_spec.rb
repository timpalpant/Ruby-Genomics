#
#  bedgraph_spec.rb
#  BioRuby
#
#  Created by Timothy Palpant on 4/8/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'bedgraph'

describe BedGraphEntry do
  context "parsed from a line" do
    BEDGRAPH_ENTRY = "chr22\t1000\t5000\t960"
    
    before do
      @test = BedGraphEntry.parse(BEDGRAPH_ENTRY)
    end
    
    it "should have chromosome chr22" do
      @test.chr.should == 'chr22'
    end
    
    it "should have start 1001" do
      @test.start.should == 1001
    end
    
    it "should have stop 5000" do
      @test.stop.should == 5000
    end
    
    it "should have value 960" do
      @test.value.should == 960
    end
    
    it "should have strand +" do
      @test.strand.should == '+'
    end
  end
end

describe BedGraphFile do
  BEDGRAPH_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bedGraph')
  
  before do
    @test = BedGraphFile.new(BEDGRAPH_FILE)
  end
  
  after do
    @test.close
  end
  
  it "should allow opening with a block" do
    BedGraphFile.open(BEDGRAPH_FILE) { |bedgraph| }
  end
  
  context "when iterating over all entries" do
    it "should have 10 entries" do
      @test.count.should == 10
    end
    
    it "should return the number of skipped entries" do
      skipped = @test.each { |entry| entry }
      skipped.should == 3
    end
    
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
      BedGraphFile.foreach(BEDGRAPH_FILE) { |entrY| count += 1 }
      count.should == 10
    end
  end
  
  context "when querying for a specific chromosome" do
    it "should have 3 entries on chrI" do
      @test.count('chrI').should == 3
      
      count = 0
      @test.each('chrI') { |entry| count += 1 }
      count.should == 3
      
      @test.chr('chrI').length.should == 3
      @test['chrI'].length.should == 3
      
      count = 0
      @test.chr('chrI') { |entry| count += 1 }
      count.should == 3
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
      @test.count('chrIV', 5, 9).should == 2
    end
    
    it "should return the correct number of entries" do      
      count = 0
      @test.each('chrIV', 5, 9) { |entry| count += 1 }
      count.should == 2
    end
    
    it "should return the correct entries" do
      @test.each('chrIV', 5, 9) do |entry|
        entry.chr.should == 'chrIV'
        entry.high.should >= 5
        entry.low.should <= 9
      end
    end
  end
  
  context "when computing statistics" do
    it "should have total = 54.2" do
      @test.total.should == 54.2
    end
    
    it "should have mean = 6.775" do
      @test.mean.should == 6.775
    end
    
    it "should have standard deviation = 4.770154609653654" do
      @test.stdev.should == 4.770154609653654
    end
  end
  
  it "should allow querying for values" do
    result = @test.query('chrIII', 10, 20)
    result.start.should == 16
    result.stop.should == 20
    result.length.should == 5
    result.coverage.should == 5
    (16..20).each { |bp| result[bp].should == 2.0 }
  end
end
