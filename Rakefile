require 'bundler/setup'
require 'rspec/core/rake_task'
require './requirements_validator'

desc 'Default: run tests'
task :default => :functional_test

desc "Check for required executables"
task :check_requirements do
  validate_requirements()
end

desc "Run unit tests"
task :unit_test => [:check_requirements, :spec]

desc "Run functional tests"
task :functional_test => [:check_requirements, :test]

desc "Run specs"
RSpec::Core::RakeTask.new(:spec) do |t|
  t.pattern = 'spec/**/*_spec.rb'
  t.rspec_opts = ['-I common', '--color', '--format doc']
end

desc "Run script specs"
RSpec::Core::RakeTask.new(:test) do |t|
  t.pattern = 'test/**/*_test.rb'
  t.rspec_opts = ['-I common', '--color', '--format doc']
end
