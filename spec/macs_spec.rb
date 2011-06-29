#
#  macs_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/22/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'macs'

describe MACSEntry do
  context "parsed from a line" do
    MACS_ENTRY = "chrI\t10\t1236\t1227\t370\t8752\t3100.00\t21.04\t4.09"

    before do
      @test = MACSEntry.parse(MACS_ENTRY)
    end
    
    it "should have chromosome chrI" do
      @test.chr.should == 'chrI'
    end
    
    it "should have start 10" do
      @test.start.should == 10
    end
    
    it "should have stop 1236" do
      @test.stop.should == 1236
    end
    
    it "should have length 1227" do
      @test.length.should == 1227
    end
    
    it "should have value 21.04" do
      @test.value.should == 21.04
      @test.fold_enrichment.should == 21.04
    end
    
    it "should have strand +" do
      @test.strand.should == '+'
    end
    
    it "should have summit 380" do
      @test.summit.should == 380
    end
    
    it "should have tags 8752" do
      @test.tags.should == 8752
    end
    
    it "should have pvalue" do
      @test.pvalue.should == 10**(3100.0 / -10)
    end
    
    it "should have fdr 4.09" do
      @test.fdr.should == 4.09
    end
  end
end

describe MACSFile do
  MACS_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.macs')
  
  before do
    @test = MACSFile.new(MACS_FILE)
  end
  
  after do
    @test.close
  end
  
  it "should allow opening with a block" do
    MACSFile.open(MACS_FILE) { |bed| }
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
      MACSFile.foreach(MACS_FILE) { |entry| count += 1 }
      count.should == 10
    end
  end
  
  context "when querying for a specific chromosome" do
    it "should have 4 entries on chrI" do
      @test.count('chrI').should == 4
      
      count = 0
      @test.each('chrI') { |entry| count += 1 }
      count.should == 4
      
      @test.chr('chrI').length.should == 4
      @test['chrI'].length.should == 4
      
      count = 0
      @test.chr('chrI') { |entry| count += 1 }
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
      @test.count('chrXI', 70_000, 72_000).should == 2
    end
    
    it "should return the correct number of entries" do      
      count = 0
      @test.each('chrXI', 70_000, 72_000) { |entry| count += 1 }
      count.should == 2
    end
    
    it "should return the correct entries" do
      @test.each('chrXI', 70_000, 72_000) do |entry|
        entry.chr.should == 'chrXI'
        entry.high.should >= 70_000
        entry.low.should <= 72_000
      end
    end
  end
end
