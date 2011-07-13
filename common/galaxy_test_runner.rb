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
require 'cheetah'
require 'utils/unix'
require 'bio/utils/ucsc'

class GalaxyTestRunner
  CURRENT_RUBY_INTERPRETER = RbConfig::CONFIG.values_at('bindir', 'ruby_install_name').join('/')
  TEST_DATA_DIR = File.expand_path(File.dirname(__FILE__) + '/../test/test-data')
  
  def initialize(galaxy_config)
    @config = galaxy_config
  end
  
  def check_compilation
    %x[ #{CURRENT_RUBY_INTERPRETER} #{@config.path}/#{@config.command.script_name} 2>&1 ]
    return $?.success?
  end
  
  def run_test(test)
    puts "Running #{test}" if ENV['DEBUG']
    
    # Generate temp files for the outputs
    tmp_outputs = Hash.new
    test.outputs.each do |name, file| 
      # Just use Tempfile to generate temp file names
      tmp_file = Tempfile.new(name)
      path = tmp_file.path
      tmp_file.close
      
      tmp_outputs[name] = path
    end
    
    result = false
    begin
      output = %x[ #{CURRENT_RUBY_INTERPRETER} #{@config.path}/#{execute_string(test, tmp_outputs)} 2>&1 ]
      
      if not $?.success?
        puts output
        raise GalaxyTestError, "Error during script execution!"
      end
      
      result = compare_output(test.outputs, tmp_outputs)
    ensure
      tmp_outputs.each { |name, file| File.delete(file) if File.exist?(file) }
    end
    
    return result
  end
  
  private
  
  def compare_output(expected, actual, stringency = 0)
    expected.each_key do |name|
      # Automatically fail if the output file does not exist
      return false if not File.exist?(actual[name])
      
      # Expand binary files to diff them
      if @config.outputs[name].format == 'bigwig'
        wig = actual[name] + '.wig'
        begin
          Bio::Utils::UCSC.bigwig_to_wig(actual[name], wig)
          diff = File.diff(TEST_DATA_DIR+'/'+expected[name], wig)
        rescue
          raise GalaxyTestError, "Error converting bigWig to Wig for diffing!"
        ensure
          File.delete(wig) if File.exist?(wig)
        end
      else
        diff = File.diff(TEST_DATA_DIR+'/'+expected[name], actual[name])
      end
      
      return false unless diff.length <= 2*stringency
    end
    
    return true
  end
  
  # Construct a execution String to run the script
  # using the parameters in the functional test
  def execute_string(test, outputs)
    # Dictionary of variables to replace in the command string
    # (inputs and outputs)
    dictionary = Hash.new
    
    test.inputs.each do |varname, value|    
      # If this param is data, look for the file in the test data directory
      if @config.inputs.include?(varname) and @config.inputs[varname].type == 'data'
        dictionary[varname] = TEST_DATA_DIR + '/' + test.inputs[varname]
      else
        dictionary[varname] = test.inputs[varname]
      end
    end
    
    outputs.each { |varname, output_file| dictionary[varname] = output_file }
    Cheetah::Template.new(@config.command.str, dictionary).to_s.split(/\s/).join(' ')
  end
end

class GalaxyTestError < StandardError
end
