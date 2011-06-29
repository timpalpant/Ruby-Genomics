#
#  nucleosome_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/22/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'nucleosome'

describe Nucleosome do
  context "instantiated in code" do
    before do
      @test = Nucleosome.new('chrI', 100, 247)
      @test.dyad = 170
      @test.dyad_stdev = 23.4
      @test.conditional_position = 0.00012
      @test.dyad_mean = 175
      @test.occupancy = 520
    end
    
    it "should correctly output to NukeCalls format" do
      @test.to_s.should == "chrI\t100\t247\t170\t23.4\t0.00012\t175\t520"
    end
  end

  context "parsed from a line" do
    NUKE_CALL_ENTRY = "chr22\t10\t30\t20\t3.2\t0.001\t21\t10"
    
    before do
      @test = Nucleosome.parse(NUKE_CALL_ENTRY)
    end
    
    it "should have chromosome chr22" do
      @test.chr.should == 'chr22'
    end
    
    it "should have start 10" do
      @test.start.should == 10
    end
    
    it "should have stop 30" do
      @test.stop.should == 30
    end
    
    it "should have dyad 20" do
      @test.dyad.should == 20
      @test.position.should == 20
    end
    
    it "should have dyad standard deviation 3.2" do
      @test.dyad_stdev.should == 3.2
    end
    
    it "should have conditional position 0.001" do
      @test.conditional_position.should == 0.001
    end
    
    it "should have dyad mean 21" do
      @test.dyad_mean.should == 21
    end
    
    it "should have value 10" do
      @test.value.should == 10
      @test.occupancy.should == 10
    end
    
    it "should have strand +" do
      @test.strand.should == '+'
    end
        
    it "should correctly output to NukeCalls format" do
      @test.to_s.should == "chr22\t10\t30\t20\t3.2\t0.001\t21\t10.0"
    end
  end
end

describe NukeCallsFile do
  NUKE_CALLS_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.nukes')
  
  before do
    @test = NukeCallsFile.new(NUKE_CALLS_FILE)
  end
  
  after do
    @test.close
  end
  
  it "should allow opening with a block" do
    NukeCallsFile.open(NUKE_CALLS_FILE) { |nukes| }
  end
  
  it "should allow opening without a block" do
    nukes = NukeCallsFile.open(NUKE_CALLS_FILE)
    nukes.close
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
      NukeCallsFile.foreach(NUKE_CALLS_FILE) { |entry| count += 1 }
      count.should == 10
    end
  end
  
  context "when querying for a specific chromosome" do
    it "should have 4 entries on chrI" do
      @test.count('chrIV').should == 4
      
      count = 0
      @test.each('chrIV') { |entry| count += 1 }
      count.should == 4
      
      @test.chr('chrIV').length.should == 4
      @test['chrIV'].length.should == 4
      
      count = 0
      @test.chr('chrIV') { |entry| count += 1 }
      count.should == 4
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
  
  context "when querying for a specific spot id" do
    it "should raise an error" do
      lambda { @test.id('MyAwesomeSpot') }.should raise_error
    end
  end
  
  context "when computing statistics" do
    it "should have total = 638" do
      @test.total.should == 638
    end
    
    it "should have mean = 63.8" do
      @test.mean.should == 63.8
    end
    
    it "should have standard deviation = 146.00397254869472" do
      @test.stdev.should == 146.00397254869472
    end
  end
  
  it "should not allow querying for values" do
    lambda { @test.query('chrIII', 10, 20) }.should raise_error
  end
end
