#
#  compact_array_spec.rb
#  ruby-genomics
#
#  Created by Timothy Palpant on 6/25/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'spec_helper'
require 'sparse_array'

describe SparseArray do
  context "without data" do
    before do
      @test = SparseArray.new
    end
    
    it "should have nil start" do
      @test.start.should be_nil
    end
    
    it "should have nil stop" do
      @test.stop.should be_nil
    end
    
    it "should have nil length" do
      @test.length.should be_nil
    end
    
    it "should have 0 coverage" do
      @test.coverage.should == 0
    end
    
    it "should have nil sum" do
      @test.sum.should be_nil
    end
    
    it "should have nil total" do
      @test.total.should be_nil
    end
    
    it "should have nil mean" do
      @test.mean.should be_nil
    end
    
    it "should have nil stdev" do
      @test.stdev.should be_nil
    end
  end

  
  context "with data" do
    before do
      @test = SparseArray.new
      @test.set(1, 1.0)
      @test.set(2, 1.0)
      @test.set(5, 2.0)
      @test.set(10, 3.0)
    end
    
    it "should allow setting data" do
      @test.set(1_000, 4.0)
      @test[1_001] = 4.0
    end
    
    it "should raise an error when attempting to set negative indices" do
      lambda { @test.set(-3, 3) }.should raise_error
      lambda { @test[-3] = 5 }.should raise_error
    end
    
    it "should allow getting data" do
      @test.get(5).should == 2.0
      @test[5].should == 2.0
    end
    
    it "should raise an error when attempting to get data with negative indices" do
      lambda { @test.get(-4) }.should raise_error
      lambda { @test[-4] }.should raise_error
    end
    
    it "should allow querying for ranges of data" do
      @test[1..2].should == [1.0, 1.0]
      @test[1...3].should == [1.0, 1.0]
      @test[1...2].should == [1.0]
      @test[2..7].should == [1.0, nil, nil, 2.0, nil, nil]
      @test[1,3].should == [1.0, 1.0, nil]
    end
    
    it "should return an empty Array when slicing with length 0" do
      @test[4,0].should == []
    end
    
    it "should raise an error when attempting to query for a range with negative indices" do
      lambda { @test[-4..-2] }.should raise_error
      lambda { @test[3..-2] }.should raise_error
    end
    
    it "should allow iterating over values" do
      count = 0
      @test.each { |value| count += 1 }
      count.should == 10
    end
    
    it "should return the indices stored" do
      @test.indices.should == [1, 2, 5, 10]
      @test.keys.should == [1, 2, 5, 10]
    end
    
    it "should return the values stored" do
      @test.values.should == [1.0, 1.0, nil, nil, 2.0, nil, nil, nil, nil, 3.0]
    end
    
    it "should have start = 1" do
      @test.start.should == 1
      @test.low.should == 1
    end
    
    it "should have stop = 10" do
      @test.stop.should == 10
      @test.high.should == 10
    end
    
    it "should dynamically update start and stop" do
      @test[0] = 4
      @test.start.should == 0
      @test.low.should == 0
      
      @test[15] = 2
      @test.stop.should == 15
      @test.high.should == 15
    end
    
    it "should have length = 10" do
      @test.length.should == 10
    end
    
    it "should have coverage = 4" do
      @test.coverage.should == 4
    end
    
    it "should have sparsity 0.4" do
      @test.sparsity.should == 0.4
    end
    
    it "should cover indices 1..10" do
      (1..10).each { |i| @test.cover?(i).should be_true }
      @test.cover?(1,10).should be_true
    end
    
    it "should not cover other indices" do
      @test.cover?(0).should be_false
      @test.cover?(-3).should be_false
      @test.cover?(11).should be_false
      @test.cover?(nil).should be_false
      
      @test.cover?(0,4).should be_false
      @test.cover?(-4,1).should be_false
      @test.cover?(0,0).should be_false
      @test.cover?(8,11).should be_false
      @test.cover?(14, 15).should be_false
    end
    
    it "should include indices 1,2,5,10" do
      [1,2,5,10].each { |i| @test.include?(i).should be_true }
      @test.include?(1,2).should be_true
    end
    
    it "should not include other indices" do
      @test.include?(4).should be_false
      @test.include?(1,10).should be_false
      @test.include?(-1).should be_false
      @test.include?(11).should be_false
    end
    
    it "should be able to add a scalar to all elements" do
      sum = @test + 4.0
      sum[1].should == 5.0
      sum[2].should == 5.0
      sum[5].should == 6.0
      sum[10].should == 7.0
      
      sum = @test + (-4)
      sum[1].should == -3.0
      sum[2].should == -3.0
      sum[5].should == -2.0
      sum[10].should == -1.0
    end
    
    it "should be able to subtract a scalar to all elements" do
      difference = @test - 4.0
      difference[1].should == -3.0
      difference[2].should == -3.0
      difference[5].should == -2.0
      difference[10].should == -1.0
    end
    
    it "should be able to multiply a scalar to all elements" do
      product = @test * 4
      product[1].should == 4.0
      product[2].should == 4.0
      product[5].should == 8.0
      product[10].should == 12.0
      
      product = @test * -4.0
      product[1].should == -4.0
      product[2].should == -4.0
      product[5].should == -8.0
      product[10].should == -12.0
    end
    
    it "should be able to divide a scalar to all elements" do
      product = @test / 4
      product[1].should == 0.25
      product[2].should == 0.25
      product[5].should == 0.5
      product[10].should == 0.75
      
      product = @test / -4.0
      product[1].should == -0.25
      product[2].should == -0.25
      product[5].should == -0.5
      product[10].should == -0.75
    end
    
    it "should be able to add another SparseArray" do
      other = SparseArray.new
      other[1] = 3
      other[2] = 4
      other[6] = 1
      
      sum = @test + other
      sum[1].should == 4
      sum[2].should == 5
      lambda { sum[5] }.should raise_error
      lambda { sum[6] }.should raise_error
      lambda { sum[10] }.should raise_error
      sum.start.should == 1
      sum.stop.should == 2
      sum.length.should == 2
      sum.coverage.should == 2
    end
    
    it "should be able to subtract another SparseArray" do
      other = SparseArray.new
      other[1] = 3
      other[2] = 4
      other[6] = 1
      
      difference = @test - other
      difference[1].should == -2
      difference[2].should == -3
      lambda { difference[5] }.should raise_error
      lambda { difference[6] }.should raise_error
      lambda { difference[10] }.should raise_error
      difference.start.should == 1
      difference.stop.should == 2
      difference.length.should == 2
      difference.coverage.should == 2
    end
    
    it "should be able to multiply another SparseArray" do
      other = SparseArray.new
      other[1] = 3
      other[2] = 4
      other[6] = 1
      
      product = @test * other
      product[1].should == 3
      product[2].should == 4
      lambda { product[5] }.should raise_error
      lambda { product[6] }.should raise_error
      lambda { product[10] }.should raise_error
      product.start.should == 1
      product.stop.should == 2
      product.length.should == 2
      product.coverage.should == 2
    end
    
    it "should be able to divide another SparseArray" do
      other = SparseArray.new
      other[1] = 3
      other[2] = 4
      other[6] = 1
      
      product = @test / other
      product[1].should be_within(1e-14).of(1.0/3)
      product[2].should be_within(1e-14).of(0.25)
      lambda { product[5] }.should raise_error
      lambda { product[6] }.should raise_error
      lambda { product[10] }.should raise_error
      product.start.should == 1
      product.stop.should == 2
      product.length.should == 2
      product.coverage.should == 2
    end
    
    it "should compute the correct sum" do
      @test.sum.should be_within(1e-14).of(7)
    end
    
    it "should compute the correct mean" do
      @test.mean.should be_within(1e-14).of(1.75)
    end
    
    it "should compute the correct standard deviation" do
      @test.stdev.should be_within(1e-14).of(0.82915619758885)
    end
    
    it "should convert to an Array" do
      @test.to_a.should == [1.0, 1.0, nil, nil, 2.0, nil, nil, nil, nil, 3.0]
    end
    
    it "should convert to a Hash" do
      h = @test.to_hash
      h.length.should == 4
      h[1].should == 1.0
      h[2].should == 1.0
      h[5].should == 2.0
      h[10].should == 3.0
    end
  end
end
