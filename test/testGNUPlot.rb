COMMON_DIR = File.expand_path(File.dirname(__FILE__) + '/../common')
$LOAD_PATH << COMMON_DIR unless $LOAD_PATH.include?(COMMON_DIR)
require 'plot/GNUPlot'

x = Array.new(10) { |i| i * Math::PI }
y = x.map { |xi| Math.cos(xi) + rand }

GNUPlot.line_graph(	:file => '/Users/timpalpant/Desktop/test.png',
                  	:x => x,
                  	:y => y,
                  	:size => [800, 600],
                  	:title => 'Line Graph',
                  	:xlabel => 'x',
                  	:ylabel => 'y' )