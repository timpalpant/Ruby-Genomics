require 'spec_helper'
require 'fixed_precision'

describe Float do
  before do
    @test = Math::PI
  end
  
  it "should output default precision with no argument" do
    @test.to_s.should == '3.141592653589793'
  end
  
  it "should output 5 digits" do
    @test.to_s(5).should == '3.1416'
  end
  
  it "should output 0 digits" do
    @test.to_s(0).should == '3'
  end
end

describe Fixnum do
  before do
    @test = 3
  end
  
  it "should output normally with no argument" do
    @test.to_s.should == '3'
  end
  
  it "should output normally with any argument" do
    @test.to_s(2).should == '3'
  end
end