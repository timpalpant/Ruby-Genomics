require 'spec_helper'

##
# Just ensure that all scripts will compile
#
# NOTE: Cryptic Ruby global variables
# $? is the Process::Status object of the last process to terminate
# $?.success? is true if the error code returned is 0, false otherwise
##

shared_examples "scripts collection" do
  it "should compile" do
    failed = Array.new
    
    @fileset.each do |f|
      # Run the script just to make sure it compiles
      %x[ ruby1.9 #{f} 2>&1 ]
    
      # Flag this set of scripts as failed unless script compilation was successful
      failed << File.basename(f) unless $?.success?
    end

    raise "Script compilation failed! (#{failed.join(',')})" if failed.length > 0
  end
end

describe "Sequencing scripts" do
  before do
    @fileset = Dir.glob(File.expand_path(File.dirname(__FILE__)) + '/../sequencing/*.rb')
  end

  include_examples "scripts collection"
end

describe "Nucleosome scripts" do
  before do
    @fileset = Dir.glob(File.expand_path(File.dirname(__FILE__)) + '/../nucleosome/*.rb')
  end

  include_examples "scripts collection"
end

describe "Data processing scripts" do
  before do
    @fileset = Dir.glob(File.expand_path(File.dirname(__FILE__)) + '/../dataprocessing/*.rb')
  end

  include_examples "scripts collection"
end

describe "Other tool wrappers" do
  before do
    @fileset = Dir.glob(File.expand_path(File.dirname(__FILE__)) + '/../other-tools/*.rb')
  end

  include_examples "scripts collection"
end

describe "UCSC tool wrappers" do
  before do
    @fileset = Dir.glob(File.expand_path(File.dirname(__FILE__)) + '/../ucsc-tools/*.rb')
  end

  include_examples "scripts collection"
end

describe "DNA structure scripts" do
  before do
    @fileset = Dir.glob(File.expand_path(File.dirname(__FILE__)) + '/../structure/*.rb')
  end

  include_examples "scripts collection"
end

describe "WigMath scripts" do
  before do
    @fileset = Dir.glob(File.expand_path(File.dirname(__FILE__)) + '/../wigmath/*.rb')
  end

  include_examples "scripts collection"
end

describe "Visualization scripts" do
  before do
    @fileset = Dir.glob(File.expand_path(File.dirname(__FILE__)) + '/../visualization/*.rb')
  end

  include_examples "scripts collection"
end

describe "Microarray scripts" do
  before do
    @fileset = Dir.glob(File.expand_path(File.dirname(__FILE__)) + '/../microarray/*.rb')
  end

  include_examples "scripts collection"
end
