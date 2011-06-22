require 'gnuplot'

##
# Expose a module of common plotting tasks built upon GNUPlot and gnuplot gem
##
module GNUPlot
  def self.line_graph(params)
    # MANDATORY parameters
    validate_params(params, :file, :x, :y)
    filename = params[:file]
    x = params[:x]
    y = params[:y]
      
    # OPTIONAL parameters (should have defaults or be unused if nil)
    title = params[:title]
    xlabel = params[:label]
    ylabel = params[:ylabel]
    font = params[:font]
    size = params[:size] || [600, 400]
    type = params[:type] || 'png'
      
    terminal = "#{type}"
    terminal += " enhanced font \"#{font}\"" unless font.nil?
    terminal += " size #{size.first}, #{size.last}"
    
    Gnuplot.open do |gp|
      Gnuplot::Plot.new( gp ) do |plot|
        plot.term terminal
        plot.output filename
      
        plot.title title unless title.nil?
        plot.xlabel xlabel unless xlabel.nil?
        plot.ylabel ylabel unless ylabel.nil?
    
        plot.data << Gnuplot::DataSet.new( [x, y] ) do |ds|
          #ds.with = "linespoints"
          ds.notitle
        end
      end
    end
  end
  
  private
  # Raise an exception if any of the parameters in *args are not specified in params
  def self.validate_params(params, *args)
    args.each do |arg|
      raise GNUPlotError, "#{arg} is a mandatory parameter!" unless params.include?(arg)
    end
  end
end

class GNUPlotError < StandardError
end