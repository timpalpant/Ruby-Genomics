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
task :test => [:check_requirements, :spec]

desc "Run specs"
RSpec::Core::RakeTask.new do |t|
    t.rspec_opts = ['-I common', '--color', '--format doc']
end