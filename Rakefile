require 'bundler/setup'
require 'rspec/core/rake_task'
require './requirements_validator'

desc 'Default: run tests'
task :default => :test

desc "Check for required executables"
task :check_requirements do
  validate_requirements()
end

desc "Run unit and functional tests"
task :test => [:check_requirements, :unit, :functional]

desc "Run only unit tests"
task :unit_tests => [:check_requirements, :unit]

desc "Run functional tests"
task :functional_tests => [:check_requirements, :functional]

desc "Run functional tests"
RSpec::Core::RakeTask.new(:functional) do |t|
  t.pattern = 'test/**/*_test.rb'
  t.rspec_opts = ['-I common', '--color', '--format doc']
end

desc "Run specs"
RSpec::Core::RakeTask.new(:unit) do |t|
  t.pattern = 'spec/**/*_spec.rb'
  t.rspec_opts = ['-I common', '--color', '--format doc']
end