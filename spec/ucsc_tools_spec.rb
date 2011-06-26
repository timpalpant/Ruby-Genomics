#
#  ucsc_tools_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'ucsc_tools'

describe UCSCTrackHeader do  
  context "when instantiated in code" do
    before do
      @test = UCSCTrackHeader.new(:type => 'bedGraph', :name => 'test header', :description => 'test description', :view_limits => '0:1')
    end
    
    it "should have type bedGraph" do
      @test.type.should == "bedGraph"
    end
    
    it "should have correct name" do
      @test.name.should == 'test header'
    end
    
    it "should have correct description" do
      @test.description.should == 'test description'
    end
    
    it "should correctly output to string" do
      @test.to_s.should == "track type=bedGraph name='test header' description='test description' viewLimits=0:1"
    end
  end
  
  context "when parsed from a track line" do
    before do
      @test = UCSCTrackHeader.parse("track type=wiggle_0 name=\"My Test Wig\" description='Test Wig for Running Specs' viewLimits=0:1 autoScale=off visibility='full'")
    end
    
    it "should have type wiggle_0" do
      @test.type.should == "wiggle_0"
    end
    
    it "should have correct name" do
      @test.name.should == 'My Test Wig'
    end
    
    it "should have correct description" do
      @test.description.should == 'Test Wig for Running Specs'
    end
    
    it "should have correct viewLimits" do
      @test.view_limits.should == '0:1'
    end
    
    it "should have correct autoScale setting" do
      @test.auto_scale.should == 'off'
    end
    
    it "should have correct visibility setting" do
      @test.visibility.should == 'full'
    end
    
    it "should correctly output to string" do
      @test.to_s.should == "track type=wiggle_0 name='My Test Wig' description='Test Wig for Running Specs' autoScale=off visibility=full viewLimits=0:1"
    end
  end
end
