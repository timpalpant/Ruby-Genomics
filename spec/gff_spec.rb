#
#  gff_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/22/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'gff'

describe GFFEntry do
  context "instantiated in code" do
    before do
      @test = GFFEntry.new('chrI', 10, 5, 'myspot', 2.0)
      @test.source = 'Nimblegen'
      @test.feature = 'ORF'
    end
    
    it "should correctly output to GFF format" do
      @test.to_gff.should == "chrI\tNimblegen\tORF\t5\t10\t2.0\t-\t.\tprobe_id=myspot;count=1"
    end
  end

  context "parsed from a line" do
    GFF_ENTRY = "chr22\tNimblegen\tORF\t1000\t5000\t960\t-\t.\tprobe_id=cloneA;count=1"
    
    before do
      @test = GFFEntry.parse(GFF_ENTRY)
    end
    
    it "should have chromosome chr22" do
      @test.chr.should == 'chr22'
      @test.seqname.should == 'chr22'
    end
    
    it "should have source Nimblegen" do
      @test.source.should == 'Nimblegen'
    end
    
    it "should have feature ORF" do
      @test.feature.should == 'ORF'
    end
    
    it "should have start 5000" do
      @test.start.should == 5000
    end
    
    it "should have stop 1000" do
      @test.stop.should == 1000
    end
    
    it "should have id/group cloneA" do
      @test.id.should == 'cloneA'
      @test.group.should == 'cloneA'
    end
    
    it "should have value 960" do
      @test.value.should == 960
      @test.score.should == 960
    end
    
    it "should have strand -" do
      @test.strand.should == '-'
    end
        
    it "should correctly output to GFF format" do
      @test.to_gff.should == "chr22\tNimblegen\tORF\t1000\t5000\t960.0\t-\t.\tprobe_id=cloneA;count=1"
    end
  end
end

describe GFFFile do
  GFF_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.gff')
  
  before do
    @test = GFFFile.new(GFF_FILE)
  end
  
  after do
    @test.close
  end
  
  it "should allow opening with a block" do
    GFFFile.open(GFF_FILE) { |gff| }
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
      GFFFile.foreach(GFF_FILE) { |entry| count += 1 }
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
    it "should return a spot with a unique id" do
      spot = @test.id('Spot2')
      spot.id.should == 'Spot2'
      spot.chr.should == 'chrI'
      spot.start.should == 101
      spot.stop.should == 95
      spot.value.should == 13.2
    end

    it "should raise an error if a specific id is not found" do
      lambda { @test.id('MyAwesomeSpot') }.should raise_error
    end
    
    it "should raise an error if multiple spots share the same id" do
      lambda { @test.id('Spot1') }.should raise_error
    end
  end
  
  context "when computing statistics" do
    it "should have total = 54.2" do
      @test.total.should == 54.2
    end
    
    it "should have mean = 6.775" do
      puts @test.mean
      @test.mean.should == 6.775
    end
    
    it "should have standard deviation = 4.770154609653654" do
      puts @test.stdev
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
