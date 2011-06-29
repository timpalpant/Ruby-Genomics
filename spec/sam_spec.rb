require 'spec_helper'
require 'sam'
require 'bio'

describe SAMEntry do  
  context "single-end crick entry" do
    # This is a Crick (reverse-complement) entry
    SINGLE_END_ENTRY = "SRR060808.16	16	chrI	24952	255	36M	*	0	0	TTATAAATCTGGTGCGACAGCTTATATTAATAAAGC	+:I4II9:CIIII6A%IIBIBIICIIIIIIIIIIII	XA:i:2	MD:Z:15T12T7	NM:i:2"
    
    before(:each) do
      @test = SAMEntry.parse(SINGLE_END_ENTRY)
    end
    
    it "is crick" do
      @test.should_not be_watson
      @test.should be_crick
    end
    
    it "is on chromosome chrI" do
      @test.chr.should == 'chrI'
      @test.rname.should == 'chrI'
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
    
    it "should have pos 24952" do
      @test.pos.should == 24952
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
    
    it "should have the correct sequence " do
      @test.seq.should == Bio::Sequence::NA.new("TTATAAATCTGGTGCGACAGCTTATATTAATAAAGC").reverse_complement
    end
    
    it "should have the correct quality score" do
      @test.qual.should == "+:I4II9:CIIII6A%IIBIBIICIIIIIIIIIIII".reverse
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
    
    it "is on chromosome chrII" do
      @test.chr.should == 'chrII'
      @test.rname.should == 'chrII'
    end
    
    it "has start 26958" do
      @test.start.should == 26958
    end
    
    it "has stop 27080" do
      @test.stop.should == 27080
    end
    
    it "should have pos 26958" do
      @test.pos.should == 26958
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
    
    it "should have the correct sequence" do
      @test.seq.should == Bio::Sequence::NA.new("ATACATAGTCTCCAGGTTGGTAAAGATGAGTCTTA")
    end
    
    it "should have the correct quality score" do
      @test.qual.should == "###################################"
    end
    
    it "should have equal length sequence and quality score" do
      @test.seq.length.should == @test.qual.length
    end
  end
end

# test.sam and test.bam in the fixtures have the same entries
# so many examples can be shared
shared_examples "sam file" do
  context "when iterating over all entries" do
    it "should correctly skip comment lines" do
      count = 0
      @test.each { |entry| count += 1 if entry.qname.start_with?('@') }
      count.should == 0
    end
    
    it "has 69 entries" do
      @test.count.should == 69
      
      count = 0
      @test.each { |entry| count += 1 }
      count.should == 69
    end
    
    it "has 43 single-end entries" do
      count = 0
      @test.each { |entry| count += 1 if entry.single? }
      count.should == 43
    end
    
    it "has 26 paired-end entries" do
      count = 0
      @test.each { |entry| count += 1 if entry.paired? }
      count.should == 26
    end
    
    it "should have 4 chromosomes" do
      @test.chromosomes.length.should == 16
    end
  end
  
  context "when iterating over reads" do   
    it "should correctly skip comment lines" do
      count = 0
      @test.each_read { |read| count += 1 if read.qname.start_with?('@') }
      count.should == 0
    end
    
    it "should have 56 reads" do
      @test.count_reads.should == 56
      
      count = 0
      @test.each_read { |read| count += 1 }
      count.should == 56
    end
  end
  
  context "when querying for a specific chromosome" do
    it "should have 3 entries on chrI" do
      @test.count('chrII').should == 3
      
      count = 0
      @test.each('chrII') { |entry| count += 1 }
      count.should == 3
      
      @test.chr('chrII').length.should == 3
      @test['chrII'].length.should == 3
      
      count = 0
      @test.chr('chrII') { |entry| count += 1 }
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
      @test.count('chrII', 582_000, 582_460).should == 1
    end
    
    it "should return the correct number of entries" do      
      count = 0
      @test.each('chrII', 582_000, 582_460) { |entry| count += 1 }
      count.should == 1
    end
    
    it "should return the correct entry" do
      @test.each('chrII', 582_000, 582_460) do |entry|
        entry.chr.should == 'chrII'
        entry.start.should == 582488
        entry.stop.should == 582453
        entry.seq.should == Bio::Sequence::NA.new("gcaagaaactgcccatcaggagatattttcgcacaa")
      end
    end
  end
end

describe SAMFile do
  TEST_SAM = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.sam')
  
  before do
    @test = SAMFile.new(TEST_SAM)
  end
  
  after do
    @test.close
  end
  
  it "should allow opening with a block" do
    SAMFile.open(TEST_SAM) { |sam| }
  end
  
  it "should allow opening without a block" do
    sam = SAMFile.open(TEST_SAM)
    sam.close
  end
  
  it "should allow static iteration over entries" do
    count = 0
    SAMFile.foreach(TEST_SAM) { |entry| count += 1 }
    count.should == 69
  end
  
  it "should allow static iteration over reads" do
    count = 0
    SAMFile.foreach_read(TEST_SAM) { |entry| count += 1 }
    count.should == 56
  end
  
  it "should return the number of skipped entries" do
    skipped = @test.each { |entry| entry }
    skipped.should == 2
  end
  
  include_examples "sam file"
end

describe BAMFile do
  TEST_BAM = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.bam')
  
  before do
    @test = BAMFile.new(TEST_BAM)
  end
  
  after do
    @test.close
  end
  
  it "should allow opening with a block" do
    BAMFile.open(TEST_BAM) { |bam| }
  end
  
  it "should allow opening without a block" do
    bam = BAMFile.open(TEST_BAM)
    bam.close
  end
  
  it "should allow static iteration over entries" do
    count = 0
    BAMFile.foreach(TEST_BAM) { |entry| count += 1 }
    count.should == 69
  end
  
  it "should allow static iteration over reads" do
    count = 0
    BAMFile.foreach_read(TEST_BAM) { |entry| count += 1 }
    count.should == 56
  end
  
  include_examples "sam file"
end
