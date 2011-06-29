#
#  spot_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'spot'

describe Spot do  
  context "with an ID" do
    before do
      @test = Spot.new("chrI", 1, 100, "MySpot")
    end
    
    it "should output to Bed format correctly" do
      @test.to_bed.should == "chrI\t0\t100\tMySpot\t.\t+"
    end
    
    it "should output to BedGraph format correctly" do
      @test.to_bedgraph.should == "chrI\t0\t100"
    end
    
    it "should output to GFF format correctly" do
      @test.to_gff.should == "chrI\tSpotArray\tfeature\t1\t100\t.\t+\t.\tprobe_id=MySpot;count=1"
    end
  end
  
  context "with an ID and a value" do
    before do
      @test = Spot.new("chrI", 1, 100, "MySpot", 10.0)
    end
    
    it "should output to Bed format correctly" do
      @test.to_bed.should == "chrI\t0\t100\tMySpot\t10.0\t+"
    end
    
    it "should output to BedGraph format correctly" do
      @test.to_bedgraph.should == "chrI\t0\t100\t10.0"
    end
    
    it "should output to GFF format correctly" do
      @test.to_gff.should == "chrI\tSpotArray\tfeature\t1\t100\t10.0\t+\t.\tprobe_id=MySpot;count=1"
    end
  end
  
  context "with a value but no ID" do
    before do
      @test = Spot.new("chrI", 55, 50, nil, 10.0)
    end
    
    it "should output to Bed format correctly" do
      @test.to_bed.should == "chrI\t49\t55\t.\t10.0\t-"
    end
    
    it "should output to BedGraph format correctly" do
      @test.to_bedgraph.should == "chrI\t49\t55\t10.0"
    end
    
    it "should output to GFF format correctly" do
      @test.to_gff.should == "chrI\tSpotArray\tfeature\t50\t55\t10.0\t-\t.\tprobe_id=no_id;count=1"
    end
  end
end
