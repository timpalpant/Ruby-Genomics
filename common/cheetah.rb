#
#  cheetah.rb
#  ruby-genomics
#  Parse Python Cheetah Templates
#  For documentation, see: http://www.cheetahtemplate.org/
#  NOTE: This is a crude implentation, only as much as needed to 
#  parse and execute Galaxy functional tests 
#
#  Created by Timothy Palpant on 7/12/11.
#  Copyright 2011 UNC. All rights reserved.
#

module Cheetah
  class Template
    def initialize(str, dictionary = nil)
      @str = str.clone
      @dictionary = dictionary.clone
    end
    
    # Remove comments
    def remove_comments
      @str.gsub!(/##.*$/, '')
    end
    
    # Parse for each loops
    def for_loop
      control_block('#for', '#end for') do |block|        
        tokens = block.split(/\s/)
        element_name = tokens[1].delete('$').delete('{').delete('}')
        array_name = tokens[3].delete('$').delete('{').delete('}')
        loop_code = tokens[4..-3].join(' ') # TODO: Is joining with spaces okay?

        # Execute the loop on the dictionary values
        result = @dictionary[array_name].map do |value|
          # Recurse, with the current element as the dictionary
          Template.new(loop_code, {element_name => value}).resolve
        end
        
        result.join(' ')
      end
    end
    
    # Parse if statements
    def if_statement
      control_block('#if', '#end if') do |block|
        # TODO: Support one-line if statements
        tokens = block[0..-8].split(/$/)
        condition = tokens[0][4..-1]
        code_block = tokens[1..-1].join("\n")
        
        # Eval the conditional statement
        if eval(Template.new(condition, @dictionary).resolve)
          Template.new(code_block, @dictionary).evaluate
        else
          String.new
        end
      end
    end
    
    # Replace variables with dictionary values
    def replace_variables
      @str.gsub!(/[$]{?\w*[\b\w}]/) do |match|
        varname = match[1..-1].delete('{').delete('}')
        
        if @dictionary.include?(varname)
          @dictionary[varname]
        else
          match
        end
      end
    end
    
    # Evaluate and parse control structures
    def evaluate
      for_loop
      if_statement
      @str
    end

    # Resolve variables with values from the dictionary
    def resolve
      replace_variables
      @str
    end
    
    # This Template string, with all dictionary elements filled in
    # and control statements parsed
    def to_s
      remove_comments
      evaluate
      resolve
    end
    
    private
    
    # Parse a control statement
    def control_block(search_start, search_stop)
      while (start = @str.index(search_start))
        stop = @str.index(search_stop, start) + search_stop.length
        @str[start..stop] = yield @str[start..stop]
      end
    end
  end
end
