require 'spec_helper'
require 'math_utils'

describe Math do
	A = 1
	B = -4
	C = 2
	
	it "should correctly return minimum value" do
		Math.min(A,B).should == B
		Math.min(B,A).should == B
		Math.min(B,C).should == B
		Math.min(C,B).should == B
		Math.min(A,C).should == A
		Math.min(C,A).should == A
	end
	
	it "should correctly return maximum value" do
		Math.max(A,B).should == A
		Math.max(B,A).should == A
		Math.max(B,C).should == C
		Math.max(C,B).should == C
		Math.max(A,C).should == C
		Math.max(C,A).should == C
	end
	
	it "should correctly return values that are equal" do
		Math.max(A,A).should == A
		Math.min(A,A).should == A
	end
	
	it "should return nil if either value is nil" do
		Math.min(A,nil).should == nil
		Math.max(A,nil).should == nil
		Math.min(nil,B).should == nil
		Math.max(nil,B).should == nil
	end
end