#
#  galaxy_config.rb
#  ruby-genomics
#  Parse Galaxy tool configuration XML files
#  For specification, see: https://bitbucket.org/galaxy/galaxy-central/wiki/ToolConfigSyntax
#  NOTE: This is not a complete implementation at the moment,
#  it is only designed for parsing and running the tests portion
#
#  Created by Timothy Palpant on 6/30/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'rexml/document'

class GalaxyConfig
  attr_accessor :path, :id, :name, :version, :description, :command, :help
  attr_reader :inputs, :outputs, :tests

  def initialize(path, id, name, description = nil, version = nil)
    @path = path
    @id = id
    @name = name
    @description = description
    @version = version
    
    @inputs = Hash.new
    @outputs = Hash.new
    @tests = Array.new
  end
  
  # Parse Galaxy config from a file
  def self.load(config_file)
    # Load the XML document
    path = File.dirname(File.expand_path(config_file))
    xml = nil
    File.open(File.expand_path(config_file)) do |f|
      begin
        xml = REXML::Document.new(f)
      rescue
        raise GalaxyConfigError, "Could not parse Galaxy config XML!"
      end
    end
    
    # Parse the tool info
    id = xml.root.attributes['id'] if xml.root.attributes.include?('id')
    name = xml.root.attributes['name'] if xml.root.attributes.include?('name')
    version = xml.root.attributes['version'] if xml.root.attributes.include?('version')
    description = xml.root.elements['description'].text if xml.root.elements['description']
    
    config = self.new(path, id, name, description, version)
    
    # Command info
    config.command = GalaxyCommand.parse(xml.root.elements['command'])
    
    # Help
    config.help = xml.root.elements['help'].text if xml.root.elements['help']
    
    # Inputs
    if xml.root.elements['inputs']
      xml.root.elements['inputs'].elements.each do |e|
        next unless e.name == 'param'
        input = GalaxyInput.parse(e)
        config.inputs[input.name] = input
      end
    end
    
    # Outputs
    if xml.root.elements['outputs']
      xml.root.elements['outputs'].elements.each do |e|
        next unless e.name == 'data'
        output = GalaxyOutput.parse(e)
        config.outputs[output.name] = output
      end
    end

    # Tests
    if xml.root.elements['tests']
      xml.root.elements['tests'].elements.each do |e|
        next unless e.name == 'test'
        config.tests << GalaxyTest.parse(e)
      end
    end
    
    return config
  end

  # Return the name and description (if it exists)
  def long_name
    if @description
      @name + ' ' + @description
    else
      @name
    end
  end
end

class GalaxyCommand
  attr_accessor :str, :interpreter
  
  def initialize(str, interpreter = nil)
    @str = str
    @interpreter = interpreter
  end
  
  def self.parse(e)
    str = e.text
    interpreter = e.attributes['interpreter']
    return self.new(str, interpreter)
  end

  def script_name
    @str.split(' ').first
  end
end

class GalaxyInput
  attr_accessor :name, :type, :label, :formats, :value
  
  def initialize(name, type, label)
    @name = name
    @type = type
    @label = label
  end
  
  def self.parse(e)
    name = e.attributes['name']
    type = e.attributes['type']
    label = e.attributes['label']
    label = e.text if label.nil?
  
    input = self.new(name, type, label)
  
    # TODO: Parse all input attributes
    case type
    when 'data'
      input.formats = e.attributes['format'].split(',') if e.attributes.include?('format')
    end

    input.value = e.attributes['value'] if e.attributes.include?('value')
  
    return input
  end
end

class GalaxyOutput
  attr_accessor :name, :format, :metadata_source
  
  def initialize(name, format)
    @name = name
    @format = format
  end
  
  def self.parse(e)
    name = e.attributes['name']
    format = e.attributes['format']
    output = self.new(name, format)
  
    # TODO: Parse all output attributes
    output.metadata_source = e.attributes['metadata_source']
  
    return output
  end
end

class GalaxyTest
  attr_accessor :inputs, :outputs
  
  def initialize(inputs = {}, outputs = {})
    @inputs = inputs
    @outputs = outputs
  end
  
  def self.parse(e)
    inputs = Hash.new
    outputs = Hash.new
  
    e.elements.each do |node|
      if node.name == 'param'
        inputs[node.attributes['name']] = node.attributes['value']
      elsif node.name == 'output'
        outputs[node.attributes['name']] = node.attributes['file']
      end
    end

    return self.new(inputs, outputs)
  end
  
  def to_s
    str = "GalaxyTest: "
    str += inputs.map { |name, value| "#{name}:#{value}" }.join(', ') + '; '
    str += outputs.map { |name, file| "#{name}:#{file}" }.join(', ')
    
    return str
  end
end

class GalaxyConfigError < StandardError
end
