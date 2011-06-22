require 'spec_helper'
require 'sam'

TEST_SAM = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.sam')

describe SAMEntry do
  context "single-end entry" do
    # This is a Crick (reverse-complement) entry
    SINGLE_END_ENTRY = "SRR060808.16	16	chrI	24952	255	36M	*	0	0	TTATAAATCTGGTGCGACAGCTTATATTAATAAAGC	+:I4II9:CIIII6A%IIBIBIICIIIIIIIIIIII	XA:i:2	MD:Z:15T12T7	NM:i:2"
    
    before(:each) do
      @test = SAMEntry.parse(SINGLE_END_ENTRY)
    end
    
    it "is crick" do
      @test.watson?.should be_false
      @test.crick?.should be_true
    end
    
    it "has length 36" do
      @test.length.should == 36
    end
    
    it "has start 24987" do
      @test.start.should == 24987
    end
    
    it "has stop 24952" do
      @test.stop.should == 24952
    end
    
    it "should be single-end" do
      @test.paired?.should be_false
      @test.unpaired?.should be_true
      @test.single?.should be_true
    end
    
    it "is mapped" do
      @test.mapped?.should be_true
      @test.unmapped?.should be_false
    end
    
    it "is a primary mapping" do
      @test.primary?.should be_true
    end
    
    it "should have equal-length sequence and quality score" do
      @test.seq.length.should == @test.qual.length
    end
  end
  
  context "watson paired-end entry" do
    PAIRED_END_ENTRY = "UNC1-RDR301647_0015:1:1:1093:13632#GCCAAT	163	chrII	26958	255	35M	=	27053	123	ATACATAGTCTCCAGGTTGGTAAAGATGAGTCTTA	###################################	XA:i:0	MD:Z:35	NM:i:0"
    
    before(:each) do
      @test = SAMEntry.parse(PAIRED_END_ENTRY)
    end
    
    it "is watson" do
      @test.watson?.should be_true
      @test.crick?.should be_false
    end
    
    it "has length 123" do
      @test.length.should == 123
    end
    
    it "has start 26958" do
      @test.start.should == 26958
    end
    
    it "has stop 27080" do
      @test.stop.should == 27080
    end
    
    it "should be paired-end" do
      @test.paired?.should be_true
      @test.unpaired?.should be_false
      @test.single?.should be_false
    end
    
    it "is mapped" do
      @test.mapped?.should be_true
      @test.unmapped?.should be_false
    end
    
    it "is a primary mapping" do
      @test.primary?.should be_true
    end
    
    it "should have equal length sequence and quality score" do
      @test.seq.length.should == @test.qual.length
    end
  end
end

describe SAMFile do
  context "iterating over entries" do
    it "should correctly skip comment lines" do
      count = 0
      SAMFile.foreach(TEST_SAM) { |entry| count += 1 if entry.qname.start_with?('@') }
      count.should == 0
    end
    
    it "has 69 entries" do
      count = 0
      SAMFile.foreach(TEST_SAM) { |entry| count += 1 }
      count.should == 69
    end
    
    it "has 43 single-end entries" do
      count = 0
      SAMFile.foreach(TEST_SAM) { |entry| count += 1 if entry.single? }
      count.should == 43
    end
    
    it "has 26 paired-end entries" do
      count = 0
      SAMFile.foreach(TEST_SAM) { |entry| count += 1 if entry.paired? }
      count.should == 26
    end
  end
  
  context "iterating over reads" do   
    it "should correctly skip comment lines" do
      count = 0
      SAMFile.foreach_read(TEST_SAM) { |read| count += 1 if read.qname.start_with?('@') }
      count.should == 0
    end
    
    it "has 56 reads" do
      count = 0
      SAMFile.foreach_read(TEST_SAM) { |read| count += 1 }
      count.should == 56
    end
  end
end
