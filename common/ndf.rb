# TODO: Load all NDF data
# Currently, just a hash between array XY-coordinates and chromosomal loci
class Ndf < Hash
	def self.load(file)
		ndf = self.new

		File.foreach(file) do |line|
			# Skip header
			next if line[0..4] == 'PROBE'

			record = line.split("\t")
			x = record[15].to_i
			y = record[16].to_i
			ndf[[x,y]] = record[12]
		end

		return ndf
	end

end