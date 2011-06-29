require 'spec_helper'
require 'genomic_interval'

describe GenomicInterval do
  CHROM = 'chrI'
  
  context "with Watson coordinates" do 
    WATSON_START = 180
    WATSON_STOP = 200
    
    before do
      @test = GenomicInterval.new(CHROM, WATSON_START, WATSON_STOP)
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
    
    it "should have center 190" do
      @test.center.should == 190
    end
    
    it "should be watson" do
      @test.should be_watson
    end
    
    it "should not be crick" do
      @test.should_not be_crick
    end
    
    it "should be on the + strand" do
      @test.strand.should == '+'
    end
    
    it "should be valid" do
      @test.should be_valid
    end
    
    it "should output to Bed format correctly" do
      @test.to_bed.should == "#{CHROM}\t#{WATSON_START-1}\t#{WATSON_STOP}\t.\t.\t+"
    end
    
    it "should output to BedGraph format correctly" do
      @test.to_bedgraph.should == "#{CHROM}\t#{WATSON_START-1}\t#{WATSON_STOP}"
    end
    
    it "should output to GFF format correctly" do
      @test.to_gff.should == "#{CHROM}\tSpotArray\tfeature\t#{WATSON_START}\t#{WATSON_STOP}\t.\t+\t.\tprobe_id=no_id;count=1"
    end
  end
  
  context "with Crick coordinates" do
    CRICK_START = 200
    CRICK_STOP = 180
    
    before do
      @test = GenomicInterval.new(CHROM, CRICK_START, CRICK_STOP)
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
    
    it "should have center 190" do
      @test.center.should == 190
    end
    
    it "should not be watson" do
      @test.should_not be_watson
    end
    
    it "should be crick" do
      @test.should be_crick
    end
    
    it "should be on the - strand" do
      @test.strand.should == '-'
    end
    
    it "should be valid" do
      @test.should be_valid
    end
    
    it "should output to Bed format correctly" do
      @test.to_bed.should == "#{CHROM}\t#{CRICK_STOP-1}\t#{CRICK_START}\t.\t.\t-"
    end
    
    it "should output to BedGraph format correctly" do
      @test.to_bedgraph.should == "#{CHROM}\t#{CRICK_STOP-1}\t#{CRICK_START}"
    end
    
    it "should output to GFF format correctly" do
      @test.to_gff.should == "#{CHROM}\tSpotArray\tfeature\t#{CRICK_STOP}\t#{CRICK_START}\t.\t-\t.\tprobe_id=no_id;count=1"
    end
  end
  
  context "with null coordinates" do
    before do
      @test = GenomicInterval.new
    end
    
    it "should have length nil" do
      @test.length.should be_nil
    end
    
    it "should have low nil" do
      @test.low.should be_nil
    end
    
    it "should have high nil" do
      @test.high.should be_nil
    end
    
    it "should have center nil" do
      @test.center.should be_nil
    end
  end
  
  context "with invalid start coordinate" do
    before do
      @test = GenomicInterval.new(CHROM, -30, 4)
    end
    
    it "should be invalid" do
      @test.should_not be_valid
    end
  end
  
  context "with invalid stop coordinate" do
    before do
      @test = GenomicInterval.new(CHROM, 9, -1)
    end
    
    it "should be invalid" do
      @test.should_not be_valid
    end
  end
  
  context "with zero start coordinate" do
    before do
      @test = GenomicInterval.new(CHROM, 0, 10)
    end
    
    it "should be invalid" do
      @test.should_not be_valid
    end
  end
  
end
