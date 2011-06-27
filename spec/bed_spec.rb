#
#  bed_spec.rb
#  BioRuby
#
#  Created by Timothy Palpant on 4/8/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'bed'

describe BedEntry do
  context "instantiated in code" do
    before do
      @test = BedEntry.new('chrI', 10, 5, 'myspot', 2.0)
    end
    
    it "should have name = id" do
      @test.name.should == 'myspot'
    end
    
    it "should correctly output to Bed-6 format" do
      @test.to_bed.should == "chrI\t4\t10\tmyspot\t2.0\t-"
    end
    
    it "should correctly output to Bed-12 format if necessary" do      
      @test.thick_start = 4
      @test.to_bed.should == "chrI\t4\t10\tmyspot\t2.0\t-\t3\t0\t0"
    end
  end

  context "parsed from a line" do
    BED_ENTRY = "chr22\t1000\t5000\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488,\t0,3512"

    before do
      @test = BedEntry.parse(BED_ENTRY)
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
    
    it "should have id/name cloneA" do
      @test.id.should == 'cloneA'
      @test.name.should == 'cloneA'
    end
    
    it "should have value 960" do
      @test.value.should == 960
    end
    
    it "should have strand +" do
      @test.strand.should == '+'
    end
    
    it "should have thickStart 1001" do
      @test.thick_start.should == 1001
    end
    
    it "should have thickEnd 5000" do
      @test.thick_end.should == 5000
    end
    
    it "should have itemRGB 0" do
      @test.item_rgb.should == '0'
    end
    
    it "should have blockCount 2" do
      @test.block_count.should == 2
    end
    
    it "should have blockSizes 567,488" do
      @test.block_sizes.should == [567, 488]
    end
    
    it "should have blockStarts 1001,4513" do
      @test.block_starts.should == [1001, 4513]
    end
    
    it "should correctly output to Bed-12 format" do
      @test.to_bed.should == "chr22\t1000\t5000\tcloneA\t960.0\t+\t1000\t5000\t0\t2\t567,488\t0,3512"
    end
  end
end

describe BedFile do
  BED_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bed')
  
  before do
    @test = BedFile.new(BED_FILE)
  end
  
  after do
    @test.close
  end
  
  it "should allow opening with a block" do
    BedFile.open(BED_FILE) { |bed| }
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
      BedFile.foreach(BED_FILE) { |entrY| count += 1 }
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
    
    it "should have mean = 5.42" do
      @test.mean.should == 5.42
    end
    
    it "should have standard deviation = 4.958870839213298" do
      @test.stdev.should == 4.958870839213298
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
  
  context "when writing to disk" do
    it "should output to Bed format"
    it "should output to BedGraph format"
    it "should output to GFF format"
  end
end