#
#  stats_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 5/31/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'stats'

describe Array, "that is empty" do
  before do
    @empty = Array.new
  end
  
  it "should have a sum of nil" do
    @empty.sum.should be_nil
  end
  
  it "should have a mean of nil" do
    @empty.mean.should be_nil
  end
  
  it "should have a median of nil" do
    @empty.median.should be_nil
  end
  
  it "should have stdev nil" do
    @empty.stdev.should be_nil
  end
  
  it "should have lower quartile nil" do
    @empty.lower_quartile.should be_nil
  end
  
  it "should have upper quartile nil" do
    @empty.upper_quartile.should be_nil
  end
end

describe Array, "of all zeros" do
  before do
    @zeros = Array.new(10, 0)
  end
  
  it "should have a sum of zero" do
    @zeros.sum.should == 0
  end
  
  it "should have a mean of 0" do
    @zeros.mean.should == 0
  end
  
  it "should have a median of 0" do
    @zeros.median.should == 0
  end
  
  it "should have a lower quartile of 0" do
    @zeros.lower_quartile.should == 0
  end
  
  it "should have an upper quartile of 0" do
    @zeros.upper_quartile.should == 0
  end
end

describe Array, "of 1..9" do
  before do
    @test = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    @test_shuffled = @test.shuffle
  end
  
  it "should have a sum of 45" do
    @test.sum.should == 45
    @test_shuffled.sum.should == 45
  end
  
  it "should have a mean of 5" do
    @test.mean.should == 5
    @test_shuffled.mean.should == 5
  end
  
  it "should have a median of 5" do
    @test.median.should == 5
    @test_shuffled.median.should == 5
  end
  
  it "should have a lower quartile of 2.5" do
    @test.lower_quartile.should == 2.5
    @test_shuffled.lower_quartile.should == 2.5
  end
  
  it "should have an upper quartile of 7.5" do
    @test.upper_quartile.should == 7.5
    @test_shuffled.upper_quartile.should == 7.5
  end
end

describe Array, "of 1..10" do
  before do
    @test = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    @test_shuffled = @test.shuffle
  end
  
  it "should have a sum of 55" do
    @test.sum.should == 55
    @test_shuffled.sum.should == 55
  end
  
  it "should have a mean of 5.5" do
    @test.mean.should == 5.5
    @test_shuffled.mean.should == 5.5
  end
  
  it "should have a median of 5.5" do
    @test.median.should == 5.5
    @test_shuffled.median.should == 5.5
  end
  
  it "should have a lower quartile of 3" do
    @test.lower_quartile.should == 3
    @test_shuffled.lower_quartile.should == 3
  end
  
  it "should have an upper quartile of 8" do
    @test.upper_quartile.should == 8
    @test_shuffled.upper_quartile.should == 8
  end
end
