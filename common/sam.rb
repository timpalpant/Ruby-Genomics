require 'bio'
require 'genomic_interval'
require 'genomic_data'

# An entry in a SAM file
# From the specification, see: http://samtools.sourceforge.net/SAM-1.4.pdf
# Also helpful: http://chagall.med.cornell.edu/NGScourse/SAM.pdf
# TODO: Better handling of empty (*) values
class SAMEntry < GenomicInterval
	attr_accessor :qname, :flag, :rname, :mapq, :cigar, :rnext, :pnext, :tlen, :seq, :qual
	
	def self.parse(line)
		record = line.chomp.split("\t")
	
		entry = self.new
		entry.qname = record[0]
		entry.flag = record[1].to_i
		entry.rname = record[2]
		entry.mapq = record[4].to_i
		entry.cigar = record[5]
		entry.rnext = record[6]
		entry.pnext = record[7].to_i
		entry.tlen = record[8].to_i
		entry.seq = Bio::Sequence::NA.new(record[9])
		entry.qual = record[10]
	
		extend = (entry.tlen == 0) ? entry.seq.length : entry.tlen.abs
	
		# According to the SAM specification, POS is the leftmost base
		# so 5' end for forward-mapping reads and
		# 3' end for reverse-complement mapping reads
		# Adjust start/stop appropriately
		if entry.watson?
			entry.start = record[3].to_i
			entry.stop = entry.start + extend - 1
		else
			entry.start = record[3].to_i + entry.seq.length - 1
			entry.stop = entry.start - extend + 1
		end

		return entry
	end
	
	def chr
		@rname
	end

	def pos
		@start
	end

	# FLAGS: See SAM specification
	def paired?
		(@flag & 0x0001) != 0
	end

	def unpaired?
		not paired?
	end

	def single?
		not paired?
	end

	def proper_pair?
		(@flag & 0x0002) != 0
	end

	def mapped?
		(@flag & 0x0004) == 0
	end

	def unmapped?
		not mapped?
	end

	def mate_mapped?
		(@flag & 0x0008) == 0
	end

	def mate_unmapped?
		not mate_mapped?
	end

	def watson?
		(@flag & 0x0010) == 0
	end

	def crick?
		not watson?
	end

	def mate_watson?
		(@flag & 0x0020) == 0
	end

	def mate_crick?
		not mate_watson?
	end

	def first?
		(@flag & 0x0040) != 0
	end

	def second?
		(@flag & 0x0080) != 0
	end

	def primary?
		(@flag & 0x0100) == 0
	end

	def failed?
		(@flag & 0x0200) != 0
	end

	def duplicate?
		(@flag & 0x0400) != 0
	end
end

# Load SAM files into memory
class SAM < GenomicData
	# NOTE: Only load small SAM files completely into memory
	# For larger files, use stream operations
	def self.load(filename)
		sam = self.new
	
		SAMFile.foreach(File.expand_path(filename)) do |entry|
			sam[entry.chr] ||= Array.new
			sam[entry.chr] << entry
		end

		return sam
	end

	# Only store GenomicInterval attributes (chr, start, stop) rather than entire SAM entry
	# to conserve memory, and only load one (forward) entry for paired-end reads
	def self.load_reads(filename)
		sam = self.new

		SAMFile.foreach_read(File.expand_path(filename)) do |read|
			sam[read.chr] ||= Array.new
			sam[read.chr] << GenomicInterval.new(read.start, read.stop)
		end

		return sam
	end
end

# Access SAM files through stream operations
# Best option for iterating over large SAM files
class SAMFile < File
	# Return each SAMEntry
	def self.foreach(filename)
		File.foreach(File.expand_path(filename)) do |line|
			# Skip comment lines
			next if line.start_with?('@')
			yield SAMEntry.parse(line)
		end
	end

	# Return each read, but only one (forward) entry for paired-end reads
	# Analogous to SAM.load_reads, but accessed in a stream
	def self.foreach_read(filename)
		File.foreach(File.expand_path(filename)) do |line|
			# Skip comment lines
			next if line.start_with?('@')
			entry = SAMEntry.parse(line)
			yield entry unless entry.paired? and entry.crick?
		end
	end
end