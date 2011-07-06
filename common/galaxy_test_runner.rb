#
#  galaxy_test_runner.rb
#  ruby-genomics
#  Run Galaxy functional tests and manage the output
#
#  Created by Timothy Palpant on 7/1/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'rbconfig'
require 'tempfile'
require 'galaxy_config'
require 'utils/unix'
require 'bio/utils/ucsc_tools'

module GalaxyTestRunner
  CURRENT_RUBY_INTERPRETER = RbConfig::CONFIG.values_at('bindir', 'ruby_install_name').join('/')
  TEST_DATA_DIR = File.expand_path(File.dirname(__FILE__) + '/../test/test-data')
  
  def self.check_compilation(config)
    %x[ #{CURRENT_RUBY_INTERPRETER} #{config.path}/#{config.command.script_name} 2>&1 ]
    return $?.success?
  end
  
  def self.run_tests(config)
    tests_passed = 0
    
    config.tests.each do |test|
      # Generate temp files for the outputs
      tmp_outputs = Hash.new
      test.outputs.each do |name, file| 
        # Just use Tempfile to generate temp file names
        tmp_file = Tempfile.new(name)
        path = tmp_file.path
        tmp_file.close
        
        tmp_outputs[name] = path
      end
      
      begin
        %x[ #{CURRENT_RUBY_INTERPRETER} #{config.path}/#{execute_string(config, test, tmp_outputs)} 2>&1 ]
        raise GalaxyTestError, "Error during script execution!" unless $?.success?
        tests_passed += 1 if compare_output(config, test.outputs, tmp_outputs)
      ensure
        tmp_outputs.each { |name, file| File.delete(file) if File.exist?(file) }
      end
    end
    
    return tests_passed
  end
  
  private
  
  def self.compare_output(config, expected, actual, stringency = 0)
    expected.each do |name, file|
      # Expand binary files to diff them
      if config.outputs[name].format == 'bigwig'
        wig = file + '.wig'
        UCSCTools.bigwig_to_wig(actual[name], wig)
        diff = File.diff(TEST_DATA_DIR+'/'+expected[name], wig)
        File.delete(wig)
      else
        diff = File.diff(TEST_DATA_DIR+'/'+expected[name], actual[name])
      end
      
      return false unless diff.length <= 2*stringency
    end
    
    return true
  end
  
  # Construct a String to run the script
  # using the parameters in the functional test
  def self.execute_string(config, test, outputs)
    config.command.str.gsub(/[$]{?\w+[\b\w}]/) do |match|
      varname = match[1..-1].delete('{').delete('}')
      
      # Replace variables with test parameters
      if test.inputs.include?(varname)
        # If this param is data, look for the file in the test data directory
        if config.inputs.include?(varname) and config.inputs[varname].type == 'data'
          TEST_DATA_DIR + '/' + test.inputs[varname]
        else
          test.inputs[varname]
        end
      elsif test.outputs.include?(varname)
        outputs[varname]
      else
        match
      end
    end
  end
end

class GalaxyTestError < StandardError
end
