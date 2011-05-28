COMMON_DIR = File.expand_path(File.dirname(__FILE__) + '/../common')
$LOAD_PATH << COMMON_DIR unless $LOAD_PATH.include?(COMMON_DIR)
require 'spot'

class Pos < Hash
	def self.load(file)
		pos = self.new

		File.foreach(file) do |line|
			# Skip header
			next if line[0..5] == 'SEQ_ID'

			record = line.split("\t")

			spot = Spot.new
			spot.id = record[2]
			spot.chr = record[1]
			spot.start = record[3].to_i
			spot.stop = spot.start + record[4].to_i - 1

			pos[spot.id] = spot
		end

		return pos
	end

end