require 'spec_helper'
require 'chromosome'

describe Chromosome do
	TEST_LENGTH = 10
	TEST_FILE = File.expand_path(File.dirname(__FILE__) + '/fixtures/test-chr.txt')
	# Should correspond to the data in test-chr.txt
	TEST_DATA = [0,3,4,9,0,6,44,3,5,7,8,9,5,6,3,1,1,2,3,4,5,6,13,15,18,22]
	
	context "with start = step = span = 1" do
		before do
			@test = Chromosome.new(TEST_LENGTH)
		end
		
		it "should have default start 1" do
			@test.start.should == 1
		end
		
		it "should have default step 1" do
			@test.step.should == 1
		end
		
		it "should have default span 1" do
			@test.span.should == 1
		end
		
		it "should have length #{TEST_LENGTH}" do
			@test.length.should == TEST_LENGTH
		end
		
		it "should have stop #{TEST_LENGTH}" do
			@test.stop.should == TEST_LENGTH
		end
		
		it "should include 1..#{TEST_LENGTH}" do
			for bp in 1..TEST_LENGTH
				@test.include?(bp).should be_true
			end
			
			@test.include?(0).should be_false
			@test.include?(-1).should be_false
			@test.include?(TEST_LENGTH+1).should be_false
		end
		
		it "should include ranges between 1..#{TEST_LENGTH}" do
			# Too high
			for bp in -5..TEST_LENGTH+5
				@test.include?(-1,bp).should be_false
				@test.include?(0,bp).should be_false
				@test.include?(@test.stop+1,bp).should be_false
				@test.include?(@test.stop+2,bp).should be_false
			end
			
			# Too low
			for bp in -5..TEST_LENGTH+5
				@test.include?(bp, -1).should be_false
				@test.include?(bp, 0).should be_false
				@test.include?(bp, @test.stop+1).should be_false
				@test.include?(bp, @test.stop+2).should be_false
			end
			
			# Included
			for start in 1..TEST_LENGTH
				for stop in 1..TEST_LENGTH
					@test.include?(start, stop).should be_true
				end
			end
		end
	end
	
	context "with start = 25" do
		before do
			@test = Chromosome.new(100, 25, 5, 4)
		end
		
		it "should have correct attributes" do
			@test.start.should == 25
			@test.step.should == 5
			@test.span.should == 4
			@test.stop.should == 124
			@test.length.should == 100
		end
	end
	
	context "loaded from file" do
		before do
			@test = Contig.load_wig(TEST_FILE, 5, 31)
		end
		
		it "should have correct attributes" do
			@test.start.should == 20
			@test.step.should == 5
			@test.span.should == 4
			@test.stop.should == 45
			@test.length.should == 26
		end
		
		it "should have correct data" do
			@test.length.should == TEST_DATA.length
			
			TEST_DATA.each_with_index do |value,i|
				value.should == @test[i]
			end
		end
		
		it "should account for start base pair" do
			@test.bases(20,25).length.should == 6
			@test.bases(20,25).to_a.should == [0,3,4,9,0,6]
			
			@test.bases(35,39).length.should == 5
			@test.bases(35,39).to_a.should == [1,1,2,3,4]
		end
		
		it "should allow Crick indexing" do
			@test.bases(25,20).length.should == 6
			@test.bases(25,20).to_a.should == [0,3,4,9,0,6].reverse
			
			@test.bases(39,35).length.should == 5
			@test.bases(39,35).to_a.should == [1,1,2,3,4].reverse
		end
		
		it "should correctly output to String" do
			@test.to_s.split("\n").length.should == @test.length
		end
	end
  
  context "with variableStep format" do
    
  end
end

describe Array do
	before do
		@test = Array.new(100) { rand }
	end
	
	it "should convert to chromosome with default parameters" do
		chr = @test.to_contig
		chr.start.should == 1
		chr.step.should == 1
		chr.span.should == 1
		chr.stop.should == 100
		chr.length.should == 100
		
		@test.each_with_index do |value,i|
			chr[i].should == @test[i]
		end
	end
	
	it "should allow parameters to be set" do
		chr = @test.to_contig(25, 5, 5)
		chr.start.should == 25
		chr.step.should == 5
		chr.span.should == 5
		chr.stop.should == 124
		chr.length.should == 100
		
		@test.each_index do |i|
			chr[i].should == @test[i]
		end
	end
end