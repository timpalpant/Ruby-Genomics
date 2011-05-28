COMMON_DIR = File.expand_path(File.dirname(__FILE__) + '/../common')
$LOAD_PATH << COMMON_DIR unless $LOAD_PATH.include?(COMMON_DIR)

# Load the SGD Features data
# Used for looking up
class SGDFeatures < Hash
	def initialize
		raise "Not yet implemented"
	end
end

# Represents an entry in the SGD Features database
class SGDFeature
	# Not yet implemented
end