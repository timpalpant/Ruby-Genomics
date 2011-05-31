# Encapsulate methods in the File class 
# that call native Unix command-line utilities
# such as wc, grep, head, tail...
class File	
  # Get the number of lines in a file using wc -l
  def self.num_lines(filename)
		if File.directory?(filename)
			# Recurse directories
			total = 0
			Dir.glob(filename + '/*').each do |f2|
				total += File.num_lines(f2)
			end
			return total
		end
		
    %x[ wc -l #{filename} ].split(' ').first.to_i
  end

	# Get the number of characters in a file using wc -c
	def self.num_chars(filename)
		%x[ wc -c #{filename} ].split(' ').first.to_i
	end

	# Get the number of words in a file using wc -c
	def self.num_words(filename)
		%x[ wc -w #{filename} ].split(' ').first.to_i
	end
  
  # Return an array of strings resulting from the output of grep -v
  # Alternatively, get all lines and then use Array#select
  def self.inverse_grep(filename, search_str)
    %x[ grep -v '#{search_str}' #{filename} ].split("\n")
  end
  
  # Return an array of strings resulting from the output of grep -n
  def self.grep_with_linenum(filename, search_str)
    %x[ grep -n '#{search_str}' #{filename} ].split("\n")
  end
  
  # Return an array of strings resulting from the output of grep
  def self.grep(filename, search_str)
    %x[ grep '#{search_str}' #{filename} ].split("\n")
  end
  
  # Get lines m..n of a file using tail and head
  def self.lines(filename, start_line, end_line = nil)
    if end_line.nil?   # Read to EOF
      %x[ tail -n+#{start_line} #{filename} ].split("\n")
    else   # Read a specific number of lines
      num_lines = end_line - start_line + 1
      %x[ tail -n+#{start_line} #{filename} 2>&1 | head -n #{num_lines} ].split("\n")
      # Seems to be much slower
      #%x[ sed -n '#{start_line},#{end_line}p; #{end_line+1}q' #{filename} ].split("\n")
    end
  end
  
  # Get the first n lines of a file using head
  def self.head(filename, num_lines)
    %x[ head -n #{num_lines} #{filename} ].split("\n")
  end
  
  # Get the bottom n lines of a file using tail
  def self.tail(filename, num_lines)
    %x[ tail -n#{num_lines} #{filename} ].split("\n")
  end

	# GZip a file
	def self.gzip(filename)
		%x[ gzip #{filename} ]
	end

	# Cross-platform way of finding an executable in the $PATH
	def self.which(cmd)
		exts = ENV['PATHEXT'] ? ENV['PATHEXT'].split(';') : ['']
		ENV['PATH'].split(File::PATH_SEPARATOR).each do |path|
			exts.each do |ext|
				exe = "#{path}/#{cmd}#{ext}"
				return exe if File.executable? exe
			end
		end

		return nil
	end
  
  # Concatenate files
  def self.cat(input_files, output_file)
    %x[ cat #{input_files.join(' ')} > #{output_file} ]
  end
end