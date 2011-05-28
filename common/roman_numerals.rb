require 'stats'

##
# Add methods to the Integer and String classes for working with Roman Numerals
##

VALUES = ["I", "V", "X", "L", "C", "D", "M"].reverse

NUMERAL = {
  "I" => 1,
  "V" => 5,
  "X" => 10,
  "L" => 50,
  "C" => 100,
  "D" => 500,
  "M" => 1000
}
NUMBERS = NUMERAL.values.sort.reverse

SHORTCUTS = {
  "CM" => "DCCCC",
  "CD" => "CCCC",
	"XC" => "LXXXX",
	"XL" => "XXXX",
	"IX" => "VIIII",
	"IV" => "IIII"
}
SHORTCUTS_REGEXP = Regexp.compile(SHORTCUTS.keys.join("|"))

LONGCUTS = SHORTCUTS.invert
LONGCUTS_REGEXP = Regexp.compile(LONGCUTS.keys.join("|"))
	
ROMAN_REGEXP = /^[MCDLXVI]+$/
ARABIC_REGEXP = /^[0-9]*$/

class Integer
  def to_roman
    roman = String.new

    mknum = []
    remainder = self
    NUMBERS.each do |num|
      tmp = remainder.divmod(num)
      mknum << tmp.first
      remainder = tmp.last
    end
	
  	#
    mknum.each_with_index do |num, i|
      roman << VALUES[i] * num
    end
	
		# Replace all shortcuts in the roman numeral
    return roman.gsub(LONGCUTS_REGEXP) { |match| LONGCUTS[match] }
  end
end

class String
  def parse_roman
    str = self.upcase
    raise("Not a well-formed Roman numeral!") unless str.roman?

  	# Convert all shorcuts to explicit expressions
    str.gsub!(SHORTCUTS_REGEXP) { |match| SHORTCUTS[match] }

		# Sum the total of all the roman characters
		str.split(//).map { |char| NUMERAL[char] }.sum
  end

	# Check if a string is a valid Roman numeral
  def roman?
  	# Liberal definition: all Roman numeral letters
		self.upcase =~ ROMAN_REGEXP
    # ALT: Stolen from O'Reilly's Perl Cookbook 6.23. Regular Expression Grabbag
    #self =~ /^M*(D?C{0,3}|C[DM])(L?X{0,3}|X[LC])(V?I{0,3}|I[VX])$/i
    # ALT: Catch all non-arabic numbers
    # not self.arabic?
  end
  
  # Check if a string is a valid Arabic numeral
  def arabic?
     self =~ ARABIC_REGEXP
  end
end