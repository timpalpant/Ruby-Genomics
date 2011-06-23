#
#  parallelizer.rb
#  BioRuby
#  Parallelize computations across chromosomes with ForkManager
#
#  Created by Timothy Palpant on 6/1/11.
#  Copyright 2011 UNC. All rights reserved.
#

require 'tmpdir'
require 'forkmanager'

module Enumerable
  @@pm = Parallel::ForkManager.new(2, {'tempdir' => Dir.tmpdir})
  
  def self.max_threads=(n)
    @@pm = Parallel::ForkManager.new(n.to_i, {'tempdir' => Dir.tmpdir})
  end
  
  def p_each
    self.each do |e|
      @@pm.start(e) and next
      yield(e)
      @@pm.finish(0)
    end

    @@pm.wait_all_children
  end
end
