require 'spec_helper'
require File.dirname(__FILE__) + '/../galaxy_config'
require File.dirname(__FILE__) + '/../galaxy_test_runner'

##
# Ensure that all scripts compile, and run functional tests, if available
#
# NOTE: Cryptic Ruby global variables
# $? is the Process::Status object of the last process to terminate
# $?.success? is true if the error code returned is 0, false otherwise
##

shared_examples "a scripts collection" do |dir|
  # Iterate over all of the Galaxy config files in this collection
  Dir.glob(dir + '/*.xml').each do |f|
    # Load the configuration
    config = GalaxyConfig.load(f)
    
    # Only test scripts and not wrappers for other scripts
    full_path = dir + '/' + config.command.script_name
    describe "#{config.command.script_name}: #{config.long_name}", :if => File.exist?(full_path) do
      subject { config }
      
      it "should compile" do
        GalaxyTestRunner.check_compilation(subject).should be_true
      end
      
      it "should pass #{config.tests.length} functional test(s)", :if => config.tests.length > 0 do
        GalaxyTestRunner.run_tests(subject).should == config.tests.length
      end
    end
  end
end

describe "Sequencing scripts" do
  SEQUENCING_SCRIPTS_DIR = File.expand_path(File.dirname(__FILE__) + '/../sequencing')
  it_behaves_like "a scripts collection", SEQUENCING_SCRIPTS_DIR
end

describe "Nucleosome scripts" do
  NUCLEOSOME_SCRIPTS_DIR = File.expand_path(File.dirname(__FILE__) + '/../nucleosome')
  it_behaves_like "a scripts collection", NUCLEOSOME_SCRIPTS_DIR
end

describe "Data processing scripts" do
  DATA_PROCESSING_SCRIPTS_DIR = File.expand_path(File.dirname(__FILE__) + '/../dataprocessing')
  it_behaves_like "a scripts collection", DATA_PROCESSING_SCRIPTS_DIR
end

describe "Other tool wrappers" do
  OTHER_TOOLS_SCRIPTS_DIR = File.expand_path(File.dirname(__FILE__) + '/../other-tools')
  it_behaves_like "a scripts collection", OTHER_TOOLS_SCRIPTS_DIR
end

describe "UCSC tool wrappers" do
  UCSC_TOOLS_SCRIPTS_DIR = File.expand_path(File.dirname(__FILE__) + '/../ucsc-tools')
  it_behaves_like "a scripts collection", UCSC_TOOLS_SCRIPTS_DIR
end

describe "DNA structure scripts" do
  DNA_STRUCTURE_SCRIPTS_DIR = File.expand_path(File.dirname(__FILE__) + '/../structure')
  it_behaves_like "a scripts collection", DNA_STRUCTURE_SCRIPTS_DIR
end

describe "WigMath scripts" do
  WIGMATH_SCRIPTS_DIR = File.expand_path(File.dirname(__FILE__) + '/../wigmath')
  it_behaves_like "a scripts collection", WIGMATH_SCRIPTS_DIR
end

describe "Visualization scripts" do
  VISUALIZATION_SCRIPTS_DIR = File.expand_path(File.dirname(__FILE__) + '/../visualization')
  it_behaves_like "a scripts collection", VISUALIZATION_SCRIPTS_DIR
end

describe "Microarray scripts" do
  MICROARRAY_SCRIPTS_DIR = File.expand_path(File.dirname(__FILE__) + '/../microarrays')
  it_behaves_like "a scripts collection", MICROARRAY_SCRIPTS_DIR
end
