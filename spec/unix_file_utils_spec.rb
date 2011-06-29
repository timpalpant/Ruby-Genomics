require 'spec_helper'
require 'unix_file_utils'
require 'fileutils'

describe File do
  # Fixtures for testing Unix file utils
  FILE_UTILS1 = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.file1')
  FILE_UTILS2 = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.file2')
  FILE_UTILS3 = File.expand_path(File.dirname(__FILE__) + '/fixtures/test.file3')
  
  it "should count lines" do
    File.num_lines(FILE_UTILS1).should == 13
    File.num_lines(FILE_UTILS2).should == 72
  end
  
  it "should count characters" do
    File.num_chars(FILE_UTILS1).should == 270
    File.num_chars(FILE_UTILS2).should == 10526
  end
  
  it "should count words" do
    File.num_words(FILE_UTILS1).should == 59
    File.num_words(FILE_UTILS2).should == 977
  end

  it "should inverse grep" do
    File.inverse_grep(FILE_UTILS1, 'chrIV').length.should == 9
    
    count = 0
    File.inverse_grep(FILE_UTILS1, 'chrIV') { |line| count += 1 }
    count.should == 9
  end
  
  it "should grep with line numbers" do
    File.grep_with_linenum(FILE_UTILS1, 'chrII').length.should == 3
    
    count = 0
    File.grep_with_linenum(FILE_UTILS1, 'chrII') { |n,line| count += 1 }
    count.should == 3
    
    line_numbers = Array.new
    File.grep_with_linenum(FILE_UTILS1, 'chrII') do |n,line| 
      line_numbers << n
      line.start_with?('chrII').should be_true
    end
    line_numbers.should == [3, 11, 13]
  end
  
  it "should grep" do
    File.grep(FILE_UTILS1, 'chrI').length.should == 3
    
    count = 0
    File.grep(FILE_UTILS1, 'chrI') { |line| count += 1 }
    count.should == 3
    
    File.grep(FILE_UTILS1, 'chrI') do |line|
      line.start_with?('chrI').should be_true
    end
  end

  it "should retrieve random lines" do
    File.lines(FILE_UTILS1, 1, 5).length.should == 5
    File.lines(FILE_UTILS1, 6, 8).length.should == 3
    File.lines(FILE_UTILS1, 2, 2).length.should == 1
    File.lines(FILE_UTILS1, 2, 2).first.should == "chrI	100	95	Spot2	13.2"
    File.lines(FILE_UTILS1, 10).length.should == 4
    
    count = 0
    File.lines(FILE_UTILS1, 6, 7) { |line| count += 1 }
    count.should == 2
  end
  
  it "should retrieve lines from the start of a file" do
    File.head(FILE_UTILS1, 5).length.should == 5
    
    count = 0
    File.head(FILE_UTILS1, 5) { |line| count += 1 }
    count.should == 5
  end
  
  it "should retrieve lines from the end of a file" do
    File.tail(FILE_UTILS1, 3).length.should == 3
    
    count = 0
    File.head(FILE_UTILS1, 3) { |line| count += 1 }
    count.should == 3
  end
  
  it "should gzip files" do
    gzipped = FILE_UTILS1 + '.gz'
    backup = FILE_UTILS1 + '.backup'
    begin
      FileUtils.copy(FILE_UTILS1, backup)
      File.gzip(FILE_UTILS1)
      File.exist?(gzipped).should be_true
    ensure
      FileUtils.move(backup, FILE_UTILS1)
      File.delete(gzipped) if File.exist?(gzipped)
    end
  end
  
  it "should find the location of an executable in the $PATH" do
    File.which('gzip').end_with?('gzip').should be_true
  end
  
  it "should return nil if an executable cannot be found in the $PATH" do
    File.which('nowaythatthisisarealexename').should be_nil
  end
  
  it "should add a newline at the end of a file if it doesn't exist" do
    temp_file = FILE_UTILS1 + '.tmp'
    begin
      File.newlinify(FILE_UTILS1, temp_file)
      File.num_lines(temp_file).should == File.num_lines(FILE_UTILS1) + 1
    ensure
      File.delete(temp_file) if File.exist?(temp_file)
    end
  end
  
  it "should not add a newline at the end of a file if it already exists" do
    temp_file = FILE_UTILS3 + '.tmp'
    begin
      File.newlinify(FILE_UTILS3, temp_file)
      File.num_lines(temp_file).should == File.num_lines(FILE_UTILS3)
    ensure
      File.delete(temp_file) if File.exist?(temp_file)
    end
  end
  
  it "should raise an error if attempting to concatenate a single file" do
    temp_file = FILE_UTILS3 + '.tmp'
    lambda { File.cat(FILE_UTILS3, temp_file) }.should raise_error
    File.exist?(temp_file).should be_false
  end
  
  it "should concatenate 2 files" do
    num_lines1 = File.num_lines(FILE_UTILS1)
    num_lines2 = File.num_lines(FILE_UTILS2)
    
    temp_file = FILE_UTILS1 + '.cat'
    begin
      File.cat([FILE_UTILS1, FILE_UTILS2], temp_file)
      File.num_lines(temp_file).should == num_lines1 + num_lines2
    ensure
      File.delete(temp_file) if File.exist?(temp_file)
    end
  end
  
  it "should concatenate 3 files" do
    files = [FILE_UTILS1, FILE_UTILS2, FILE_UTILS3]
    num_lines = files.map { |f| File.num_lines(f) }
    total_num_lines = num_lines[0] + num_lines[1] + num_lines[2]
    
    temp_file = FILE_UTILS1 + '.cat'
    begin
      File.cat(files, temp_file)
      File.num_lines(temp_file).should == total_num_lines
    ensure
      File.delete(temp_file) if File.exist?(temp_file)
    end
  end
  
  it "should sort files" do
    temp_file = FILE_UTILS1 + '.sorted'
    begin
      File.sort(FILE_UTILS1, temp_file, '-k1,1 -k2,2n')
      sorted_lines = File.readlines(temp_file)
      expected = File.readlines(FILE_UTILS3)
      sorted_lines.should == expected
    ensure
      File.delete(temp_file) if File.exist?(temp_file)
    end
  end
end
