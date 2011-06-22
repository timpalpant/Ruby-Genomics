require 'spec_helper'
require 'genomic_interval'

describe GenomicInterval do
  
  context "with Watson coordinates" do 
    WATSON_START = 180
    WATSON_STOP = 200
    
    before do
      @test = GenomicInterval.new
      @test.start = WATSON_START
      @test.stop = WATSON_STOP
    end
    
    it "should have length #{WATSON_STOP-WATSON_START+1}" do
      @test.length.should == (WATSON_STOP-WATSON_START+1)
    end
    
    it "should include #{WATSON_START}..#{WATSON_STOP}" do
      for bp in WATSON_START..WATSON_STOP
        @test.include?(bp).should be_true
      end
      
      # Test outside the window
      @test.include?(WATSON_START-1).should be_false
      @test.include?(WATSON_STOP+1).should be_false
      
      # Test odd numbers
      @test.include?(0).should be_false
      @test.include?(-1).should be_false
    end
    
    it "should have low = #{WATSON_START}" do
      @test.low.should == WATSON_START
    end
    
    it "should have high #{WATSON_STOP}" do
      @test.high.should == WATSON_STOP
    end
    
    it "should be watson" do
      @test.watson?.should be_true
    end
    
    it "should not be crick" do
      @test.crick?.should be_false
    end
    
    it "should be invalid" do
      @test.valid?.should be_true
    end
  end
  
  context "with Crick coordinates" do
    CRICK_START = 200
    CRICK_STOP = 180
    
    before do
      @test = GenomicInterval.new
      @test.start = CRICK_START
      @test.stop = CRICK_STOP
    end
    
    it "should have length #{CRICK_START-CRICK_STOP+1}" do
      @test.length.should == (CRICK_START-CRICK_STOP+1)
    end
    
    it "should include #{CRICK_STOP}..#{CRICK_START}" do
      for bp in CRICK_START..CRICK_STOP
        @test.include?(bp).should be_true
      end
      
      # Test outside the window
      @test.include?(CRICK_START+1).should be_false
      @test.include?(CRICK_STOP-1).should be_false
      
      # Test odd numbers
      @test.include?(0).should be_false
      @test.include?(-1).should be_false
    end
    
    it "should have low = #{CRICK_STOP}" do
      @test.low.should == CRICK_STOP
    end
    
    it "should have high #{CRICK_START}" do
      @test.high.should == CRICK_START
    end
    
    it "should not be watson" do
      @test.watson?.should be_false
    end
    
    it "should be crick" do
      @test.crick?.should be_true
    end
    
    it "should be valid" do
      @test.valid?.should be_true
    end
  end
  
  context "with null coordinates" do
    before do
      @test = GenomicInterval.new
    end
    
    it "should have length nil" do
      @test.length.should == nil
    end
    
    it "should have low nil" do
      @test.low.should == nil
    end
    
    it "should have high nil" do
      @test.high.should == nil
    end
  end
  
  context "with invalid start coordinate" do
    before do
      @test = GenomicInterval.new
      @test.start = -30
      @test.stop = 4
    end
    
    it "should be invalid" do
      @test.valid?.should be_false
    end
  end
  
  context "with invalid stop coordinate" do
    before do
      @test = GenomicInterval.new
      @test.start = 9
      @test.stop = -1
    end
    
    it "should be invalid" do
      @test.valid?.should be_false
    end
  end
  
  context "with zero start coordinate" do
    before do
      @test = GenomicInterval.new
      @test.start = 0
      @test.stop = 10
    end
    
    it "should be invalid" do
      @test.valid?.should be_false
    end
  end
  
end