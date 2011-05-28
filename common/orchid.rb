require 'perl'


##
# Do Orchid computations on sequences and FASTA files
##
module Orchid
  ORCHID_SEQUENCE_SCRIPT = 'orchid.pl'
  ORCHID_FASTA_SCRIPT = 'orchid_fasta.pl'
  ORCHID_FASTA_TEMP = 'orchid_temp.fa'
  
  # Run Orchid on a Bio::Sequence::NA
  def self.sequence(seq)
    # Use the command-line Orchid script if the sequence is < 100,000 bp
    if seq.length < 100_000
      values = Perl.run("#{ORCHID_SEQUENCE_SCRIPT} -s #{seq}").split("\n").map { |line| line.split("\t").last.to_f }
    # Otherwise, dump the sequence to a FASTA file and use the FASTA Orchid script
    else
      # Dump sequence to a temporary FASTA file
      File.open(ORCHID_FASTA_TEMP, 'w') { |f| f.puts seq.to_fasta('temp') }
      
      # Run orchid_fasta.pl on the temporary FASTA file
      values = self.fasta(ORCHID_FASTA_TEMP).split(',')[1..-1].map { |value| value.to_f }
        
      # Remove the temporary FASTA file
      File.delete(ORCHID_FASTA_TEMP)
    end
    
    raise "Number of Orchid values does not match length of sequence! (#{values.length} vs #{self.length})" if values.length != seq.length
    return values
  end
  
  # Run Orchid on a Fasta file and return an array of the output lines
  def self.fasta(filename)
    Perl.run("#{ORCHID_FASTA_SCRIPT} -f #{filename}").split("\n")
  end
end