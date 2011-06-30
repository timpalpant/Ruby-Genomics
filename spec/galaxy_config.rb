require 'rexml/document'

##
# Parse Galaxy config files (read-only)
##
class GalaxyConfig
  attr_accessor :id, :name, :version, :description, :interpreter, :command, :help
  attr_reader :inputs, :outputs, :tests

  def initialize(id, name, description = nil, version = nil)
    @id = id
    @name = name
    @description = description
    @version = version
    
    @inputs = Array.new
    @outputs = Array.new
    @tests = Array.new
  end
  
  # Parse Galaxy config
  def self.parse(io)
    # Load the XML document
    xml = REXML::Document.new(io)
    
    # Parse the tool info
    id = xml.root.attributes['id'] if xml.root.attributes.include?('id')
    name = xml.root.attributes['name'] if xml.root.attributes.include?('name')
    version = xml.root.attributes['version'] if xml.root.attributes.include?('version')
    description = xml.root.elements['description'].text if xml.root.elements['description']
    
    config = self.new(id, name, description, version)
    
    # Command info
    config.command = xml.root.elements['command'].text
    config.interpreter = xml.root.elements['command'].attributes['interpreter'] if xml.root.elements['command'].attributes.include?('interpreter')
    
    # Help
    config.help = xml.root.elements['help'].text if xml.root.elements['help']
    
    # Inputs
    xml.root.elements['inputs'].elements.each do |e|
      next unless e.name == 'param'
      config.inputs << GalaxyInput.parse(e)
    end
    
    # Outputs
    @outputs = Array.new
    xml.root.elements['outputs'].elements.each do |e|
      next unless e.name == 'data'
      config.outputs << GalaxyOutput.parse(e)
    end
    
    # Tests
    @tests = Array.new
    xml.root.elements['tests'].elements.each do |e|
      next unless e.name == 'test'
      config.tests << GalaxyTest.parse(e)
    end
    
    return config
  end
  
  # Load a Galaxy config file
  def self.load(filename)
    File.open(File.expand_path(filename)) do |f|
      config = self.parse(f)
    end
    
    return config
  end
end

class GalaxyInput
  def self.parse(e)
    
  end
end

class GalaxyOutput
  def self.parse(e)
  
  end
end

class GalaxyTest
  def self.parse(e)
  
  end
end

class GalaxyConfigError < StandardError
end
