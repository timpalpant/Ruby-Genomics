require 'spec_helper'
require 'roman_numerals'

describe Integer, "#to_roman" do
  it "should convert Roman digits correctly" do
    1.to_roman.should == 'I'
    5.to_roman.should == 'V'
    10.to_roman.should == 'X'
    50.to_roman.should == 'L'
    100.to_roman.should == 'C'
    500.to_roman.should == 'D'
    1000.to_roman.should == 'M'
  end
  
  it "should convert to proper Roman shortcuts" do
    4.to_roman.should == 'IV'
    9.to_roman.should == 'IX'
    14.to_roman.should == 'XIV'
    19.to_roman.should == 'XIX'
    49.to_roman.should == 'XLIX'
    99.to_roman.should == 'XCIX'
  end
end

describe String, "#parse_roman" do
  it "should parse Roman digits correctly" do
    'I'.parse_roman.should == 1
    'V'.parse_roman.should == 5
    'X'.parse_roman.should == 10
    'L'.parse_roman.should == 50
    'C'.parse_roman.should == 100
    'D'.parse_roman.should == 500
    'M'.parse_roman.should == 1000
  end
  
  it "should parse Roman shortcuts correctly" do
    'IV'.parse_roman.should == 4
    'IX'.parse_roman.should == 9
    'XIV'.parse_roman.should == 14
    'XIX'.parse_roman.should == 19
    'XLIX'.parse_roman.should == 49
    'XCIX'.parse_roman.should == 99
  end
  
  it "should parse improper Roman numerals correctly" do
    'IIII'.parse_roman.should == 4
    'IIIIIIIIII'.parse_roman.should == 10
    'VIIII'.parse_roman.should == 9
    'XVIIII'.parse_roman.should == 19
  end
  
  it "should raise an error if a Roman numeral is not parsable"
end

describe String, "#roman?" do
  it "should return true for Roman Numerals" do
    'I'.roman?.should be_true
    'V'.roman?.should be_true
    'X'.roman?.should be_true
    'L'.roman?.should be_true
    'C'.roman?.should be_true
    'D'.roman?.should be_true
    'M'.roman?.should be_true
  end
  
  it "should return true for compound Roman Numeral Strings"
  
  it "should return false for all other Strings"
end

describe String, "#arabic?" do
  it "should return true for all Arabic integers"
  
  it "should return false for all other Strings"
end
