PERL_EXE = 'perl'
SAMTOOLS_EXE = 'samtools'
TABIX_EXE = %w[ tabix 
                bgzip ]
UNIX_EXE = %w[ head 
               tail 
               grep 
               wc 
               sort 
               cat ]
UCSC_EXE = %w[ bedGraphToBigWig 
               bedSort 
               bedToBigBed 
               bigBedInfo 
               bigBedSummary 
               bigBedToBed 
               bigWigInfo 
               bigWigSummary 
               bigWigToBedGraph 
               bigWigToWig 
               wigCorrelate 
               wigToBigWig ]

# Cross-platform way of finding an executable in the $PATH.
# which('ruby') => /usr/bin/ruby
# which('notanexe') => nil
def which(cmd)
  exts = ENV['PATHEXT'] ? ENV['PATHEXT'].split(';') : ['']

  ENV['PATH'].split(File::PATH_SEPARATOR).each do |path|
    exts.each do |ext|
      exe = "#{path}/#{cmd}#{ext}"
      return exe if File.executable?(exe)
    end
  end

  return nil
end

def check_exe(name, exe_set)
  error_msg = "\n\n-------------------------------------------------------------------\n"
  error_msg += "Cannot find required executable(s): #{name}"
  valid = true

  if exe_set.is_a?(String)
    if which(exe_set).nil?
      valid = false
      error_msg += " (#{exe_set})\n"
    end
  elsif exe_set.is_a?(Array)
    missing = exe_set.select { |exe| which(exe).nil? }
    if missing.length > 0
      valid = false
      missing.each { |exe| error_msg += "\n\t#{exe}" }
    end
    error_msg += "\n"
  else
    raise "Do not know how to check exe specified by type #{exe_set.class}"
  end

  error_msg += "in the $PATH"
  error_msg += "\n-------------------------------------------------------------------\n\n\n"

  if not valid
    puts error_msg
    exit(-1)
  end
end

def validate_requirements
  check_exe('Perl', PERL_EXE)
  check_exe('SAMTools', SAMTOOLS_EXE)
  check_exe('Tabix', TABIX_EXE)
  check_exe('Unix Utilities', UNIX_EXE)
  check_exe('UCSC Binaries', UCSC_EXE)
end
