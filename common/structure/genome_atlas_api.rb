require 'rubygems'
require 'bio'

# Dynamically generated API for the GenomeAtlas SOAP web service
class GenomeAtlasAPI < Bio::SOAPWSDL
	def initialize
    @wsdl = 'http://www.cbs.dtu.dk/ws/GenomeAtlas/GenomeAtlas_3_0_ws2.wsdl'
    create_driver
  end
end