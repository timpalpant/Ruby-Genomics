require 'spec_helper'
require 'sam'

describe SAMEntry do
  context "single-end entry" do
    # This is a Crick (reverse-complement) entry
    SINGLE_END_ENTRY = "SRR060808.16	16	chrI	24952	255	36M	*	0	0	TTATAAATCTGGTGCGACAGCTTATATTAATAAAGC	+:I4II9:CIIII6A%IIBIBIICIIIIIIIIIIII	XA:i:2	MD:Z:15T12T7	NM:i:2"
    
    before(:each) do
      @test = SAMEntry.parse(SINGLE_END_ENTRY)
    end
    
    it "is crick" do
      @test.should_not be_watson
      @test.should be_crick
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
      @test.should_not be_paired
      @test.should be_unpaired
      @test.should be_single
    end
    
    it "is mapped" do
      @test.should be_mapped
      @test.should_not be_unmapped
    end
    
    it "is a primary mapping" do
      @test.should be_primary
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
      @test.should be_watson
      @test.should_not be_crick
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
      @test.should be_paired
      @test.should_not be_unpaired
      @test.should_not be_single
    end
    
    it "is mapped" do
      @test.should be_mapped
      @test.should_not be_unmapped
    end
    
    it "is a primary mapping" do
      @test.should be_primary
    end
    
    it "should have equal length sequence and quality score" do
      @test.seq.length.should == @test.qual.length
    end
  end
end

shared_examples "sam file" do
  
end

describe SAMFile do
  TEST_SAM = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.sam')
  
  before do
    @test_file = TEST_SAM
  end
  
  context "iterating over entries" do
    it "should correctly skip comment lines" do
      count = 0
      SAMFile.foreach(@test_file) { |entry| count += 1 if entry.qname.start_with?('@') }
      count.should == 0
    end
    
    it "has 69 entries" do
      count = 0
      SAMFile.foreach(@test_file) { |entry| count += 1 }
      count.should == 69
    end
    
    it "has 43 single-end entries" do
      count = 0
      SAMFile.foreach(@test_file) { |entry| count += 1 if entry.single? }
      count.should == 43
    end
    
    it "has 26 paired-end entries" do
      count = 0
      SAMFile.foreach(@test_file) { |entry| count += 1 if entry.paired? }
      count.should == 26
    end
  end
  
  context "iterating over reads" do   
    it "should correctly skip comment lines" do
      count = 0
      SAMFile.foreach_read(@test_file) { |read| count += 1 if read.qname.start_with?('@') }
      count.should == 0
    end
    
    it "has 56 reads" do
      count = 0
      SAMFile.foreach_read(@test_file) { |read| count += 1 }
      count.should == 56
    end
  end
end

describe BAMFile do
  TEST_BAM = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bam')
  
  before do
    @test_file = TEST_BAM
  end
  
  context "iterating over entries" do
    it "should correctly skip comment lines" do
      count = 0
      BAMFile.foreach(@test_file) { |entry| count += 1 if entry.qname.start_with?('@') }
      count.should == 0
    end
    
    it "has 69 entries" do
      count = 0
      BAMFile.foreach(@test_file) { |entry| count += 1 }
      count.should == 69
    end
    
    it "has 43 single-end entries" do
      count = 0
      BAMFile.foreach(@test_file) { |entry| count += 1 if entry.single? }
      count.should == 43
    end
    
    it "has 26 paired-end entries" do
      count = 0
      BAMFile.foreach(@test_file) { |entry| count += 1 if entry.paired? }
      count.should == 26
    end
  end
  
  context "iterating over reads" do   
    it "should correctly skip comment lines" do
      count = 0
      BAMFile.foreach_read(@test_file) { |read| count += 1 if read.qname.start_with?('@') }
      count.should == 0
    end
    
    it "has 56 reads" do
      count = 0
      BAMFile.foreach_read(@test_file) { |read| count += 1 }
      count.should == 56
    end
  end
end
