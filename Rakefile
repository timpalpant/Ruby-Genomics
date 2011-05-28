require 'rubygems'
require 'rspec/core/rake_task'

desc 'Default: run specs.'
task :default => :spec

desc "Run specs"
RSpec::Core::RakeTask.new do |t|
    t.rspec_opts = ['-I common', '--color', '--format doc']

end
