require 'genomic_data'
require 'genomic_interval'

##
# Bowtie Output Files
##
class BowtieHits < GenomicData
	# Load single-end hits file
	def self.load_single_end(filename, extend)
		reads = self.new
		
		num_reads = File.num_lines(filename)
		
		File.foreach(filename) do |line|
			next if line.start_with?('@')
			record = line.split("\t")
			
			# Get the read location
			strand = record[1]
			chr = record[2]
			start = record[3].to_i
			read = record[4]
			length = if extend.nil? or extend < 1
				read.length
			else
				extend
			end
			# Account for directionality of strand
			start = start + read.length - length if strand == '-'
			# Extend read by specified length
			stop = start + length

			# Ensure that genomic coordinates are valid
			# baseAlignCounts will do this for us
			start = 1 if start < 1

			# Add to Array
			reads[chr] ||= Array.new
			reads[chr] << GenomicInterval.new(start,stop)
		end

		return reads
	end

	# Load paired-end hits file
	def self.load_paired_end(filename)
		reads = self.new

		num_lines = File.num_lines(filename)
		num_reads = num_lines / 2

		File.open(filename) do |f|
			line_num = 0
			while line_num < num_lines
				read1 = f.gets.chomp
				read2 = f.gets.chomp
				entry1 = read1.split("\t")
				entry2 = read2.split("\t")
				
				# Ensure that the chromosomes are the same
				chr1 = entry1[2]
				chr2 = entry2[2]
				raise "Chromosomes of paired reads are mismatched!" if chr1 != chr2
				
				# Get the read start and stop locations
				start = entry1[3].to_i
				# Stop location is offset by the length of read2 
				# because of the way that Bowtie reports loci
				stop = entry2[3].to_i + entry2[4].length
				
				# Add to Array
				reads[chr] ||= Array.new
				reads[chr] << GenomicInterval.new(start,stop)
				puts line_num if line_num % 1_000_000 == 0
				line_num += 2
			end
		end

		return reads
	end
end