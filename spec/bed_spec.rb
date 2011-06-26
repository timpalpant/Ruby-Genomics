#
#  bed_spec.rb
#  BioRuby
#
#  Created by Timothy Palpant on 4/8/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'spot_file_examples'
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
    TEST_ENTRY = "chr22\t1000\t5000\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488,\t0,3512"

    before do
      @test = BedEntry.parse(TEST_ENTRY)
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
  TEST_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bed')
  
  before do
    @test = BedFile.new(TEST_FILE)
  end
  
  after do
    @test.close
  end
  
  it "should have 10 entries" do
    @test.count.should == 10
  end
  
  it "should return the number of skipped entries"
  
  it "should have 4 chromosomes" do
    @test.chromosomes
    @test.chromosomes.length.should == 4
  end
  
  it "should iterate over all the entries" do
    count = 0
    @test.each { |entry| count += 1 }
    count.should == 10
  end
  
  it "should have 3 entries on chrI" do
    count = 0
    @test.each('chrI') { |entry| count += 1 }
    count.should == 3
    
    @test.chr('chrI').length.should == 3
    
    count = 0
    @test.chr('chrI') { |entry| count += 1 }
    count.should == 3
  end
  
  it "should have 0 entries on chr8" do
    count = 0
    @test.each('chr8') { |entry| count += 1 }
    count.should == 0
    
    @test.chr('chr8').length.should == 0
    
    count = 0
    @test.chr('chr8') { |entry| count += 1 }
    count.should == 0
  end
end