require 'fileutils'

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
    
    # Add 1 because wc counts line breaks?
    %x[ wc -l #{filename} ].split(' ').first.to_i + 1
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
  # Takes an optional block
  def self.inverse_grep(filename, search_str, &block)
    if block
      IO.popen("grep -v -w '#{search_str}' #{filename}") do |output|
        output.each { |line| yield line }
      end
    else
      return %x[ grep -v '#{search_str}' #{filename} ].split("\n")
    end
  end
  
  # Return an array of strings resulting from the output of grep -n
  def self.grep_with_linenum(filename, search_str, &block)
    if block
      IO.popen("grep -n -w '#{search_str}' #{filename}") do |output|
        output.each do |line|
          entry = line.split(':')
          line_num = entry.first.to_i
          content = entry.last
          yield(line_num, content)
        end
      end
    else
      return %x[ grep -n -w '#{search_str}' #{filename} ].split("\n").map do |line|
        entry = line.split(':')
        line_num = entry.first.to_i
        content = entry.last
        [line_num, content]
      end
    end
  end
  
  # Return an array of strings resulting from the output of grep
  def self.grep(filename, search_str, &block)
    if block
      IO.popen("grep -w '#{search_str}' #{filename}") do |output|
        output.each { |line| yield line }
      end
    else
      return %x[ grep -w '#{search_str}' #{filename} ].split("\n")
    end
  end
  
  # Get lines m..n of a file using tail and head
  def self.lines(filename, start_line, end_line = nil, &block)
    raise "Cannot get lines < 1 from file!" if start_line < 1

    if end_line.nil?  # Read to EOF
      if block
        IO.popen("tail -n+#{start_line} #{filename}") do |output|
          output.each { |line| yield line }
        end
      else
        return %x[ tail -n+#{start_line} #{filename} ].split("\n")
      end
    else  # Read a specific number of lines
      num_lines = end_line - start_line + 1
      if block
        IO.popen("tail -n+#{start_line} #{filename} 2>&1 | head -n #{num_lines}") do |output|
          output.each { |line| yield line }
        end
      else
        return %x[ tail -n+#{start_line} #{filename} 2>&1 | head -n #{num_lines} ].split("\n")
        # Seems to be much slower
        #%x[ sed -n '#{start_line},#{end_line}p; #{end_line+1}q' #{filename} ].split("\n")
      end
    end
  end
  
  # Get the first n lines of a file using head
  def self.head(filename, num_lines, &block)
    if block
      IO.popen("head -n #{num_lines} #{filename}") do |output|
        output.each { |line| yield line }
      end
    else
      return %x[ head -n #{num_lines} #{filename} ].split("\n")
    end
  end
  
  # Get the bottom n lines of a file using tail
  def self.tail(filename, num_lines, &block)
    if block
      IO.popen("tail -n#{num_lines} #{filename}") do |output|
        output.each { |line| yield line }
      end
    else
      return %x[ tail -n #{num_lines} #{filename} ].split("\n")
    end
  end

  # GZip a file
  def self.gzip(filename)
    %x[ gzip #{filename} ]
  end

  # Find the location of an executable in the $PATH
  def self.which(cmd)
    %x[ which #{cmd} ].chomp
  end

  # Add a newline to the end of a file only if there isn't one already
  def self.newlinify(input_file, output_file)
    FileUtils.copy(input_file, output_file)
    if not %x[ tail -n 1 #{output_file} ].end_with?("\n")
      File.open(output_file, 'a') { |f| f.write "\n" }
    end
  end
  
  # Concatenate files
  def self.cat(input_files, output_file)
    raise "Less than 2 input files passed to cat!" if input_files.length < 2

    temp_files = input_files[0..-2].map { |f| f+'.tmp' }
    begin
      # Add newlines at the end of files 0..-2 if they don't have them
      input_files[0..-2].each { |f| File.newlinify(f, f+'.tmp') }

      # Cat all of the files together
      %x[ cat #{temp_files.join(' ')} #{input_files.last} > #{output_file} ]
    ensure
      temp_files.each { |f| File.delete(f) if File.exist?(f) }
    end
  end
  
  # Sort files
  def self.sort(input_file, output_file, options)
    %x[ sort #{options} -o #{output_file} #{input_file} ]
  end
end
