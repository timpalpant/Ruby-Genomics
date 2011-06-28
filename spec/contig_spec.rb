require 'spec_helper'
require 'contig'

describe Contig do
  before do
    @test = Contig.new('chrI')
    @test.set(100, 1.5)
    @test.set(101, 1.5)
    (105..111).each { |bp| @test.set(bp, 2.1234321) }
  end
  
  context "when no resolution is passed" do
    it "should output to variableStep format with resolution = 1" do
      @test.to_variable_step.should == "variableStep chrom=chrI span=1\n100\t1.5\n101\t1.5\n105\t2.1234\n106\t2.1234\n107\t2.1234\n108\t2.1234\n109\t2.1234\n110\t2.1234\n111\t2.1234"
    end
    
    it "should output to fixedStep format with resolution = 1" do
      @test.to_fixed_step.should == "fixedStep chrom=chrI start=100 step=1 span=1\n1.5\n1.5\nNaN\nNaN\nNaN\n2.1234\n2.1234\n2.1234\n2.1234\n2.1234\n2.1234\n2.1234"
    end
  end
  
  context "when a resolution of 2 is passed" do
    it "should output to variableStep format with resolution = 2" do
      @test.to_variable_step(2).should == "variableStep chrom=chrI span=2\n100\t1.5\n106\t2.1234\n108\t2.1234\n110\t2.1234"
    end
    
    it "should output to fixedStep format with resolution = 2" do
      @test.to_fixed_step(2).should == "fixedStep chrom=chrI start=100 step=2 span=2\n1.5\nNaN\nNaN\n2.1234\n2.1234\n2.1234"
    end
  end
end

describe Array do
  before do
    @test = [0,3,4,9,0,6,44,3,5,7,8,9,5,6,3,1,1,2,3,4,5,6,13,15,18,22]
  end
  
  it "should convert to Contig with default parameters: chr=unknown, start=1" do
    contig = @test.to_contig
    contig.chr.should == 'unknown'
    contig.start.should == 1
    contig.stop.should == @test.length
    contig.length.should == @test.length
    
    @test.each_with_index do |value,i|
      contig[i+1].should == @test[i]
    end
  end
  
  it "should allow parameters to be set" do
    contig = @test.to_contig('chrI', 25, 5, 5)
    contig.chr.should == 'chrI'
    contig.start.should == 25
    contig.stop.should == 25 + 5*(@test.length-1)+5-1
    contig.length.should == 5*(@test.length-1)+5
  end
end
