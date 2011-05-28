require 'tree'

#Just so I can raise errors with a cool names
class ValenceShellError < StandardError
end

class ElementError < StandardError
end

#Hashes for names, weights, etc., used in chemistry
module Chemistry
  ATOM_NAMES = {"H" => "hydrogen", "He" => "helium", "B" => "boron", "C" => "carbon", "N" => "nitrogen", "O" => "oxygen", "F" => "fluorine", "Ne" => "neon", "Si" => "silicon", "P" => "phosphorous", "S" => "sulfur", "Cl" => "chlorine", "Ar" => "argon", "As" => "arsenic", "Se" => "selenium", "Br" => "bromine", "Kr" => "krypton", "Te" => "tellurium", "I" => "iodine", "Xe" => "xenon", "At" => "astatine", "Rn" => "radon", "Li" => "lithium", "Na" => "sodium", "K" => "potassium", "Rb" => "rubidium", "Cs" => "cesium", "Be" => "beryllium", "Mg" => "magnesium", "Ca" => "calcium", "Sr" => "strontium", "Ba" => "barium"}
  ATOM_WEIGHTS = {"H" => 1.00794, "He" => 4.003, "B" => 10.811, "C" => 12.0107, "N" => 14.00674, "O" => 15.9994, "F" => 18.9984032, "Ne" => 20.1797, "Si" => 28.0855, "P" => 30.973761, "S" => 32.066, "Cl" => 35.4527, "Ar" => 39.948, "As" => 74.92160, "Se" => 78.96, "Br" => 79.904, "Kr" => 83.80, "Te" => 127.60, "I" => 126.90447, "Xe" => 131.29, "At" => 210.0, "Rn" => 222.0, "Li" => 6.941, "Na" => 22.989770, "K" => 39.0983, "Rb" => 85.4678, "Cs" => 132.90545, "Be" => 9.012182, "Mg" => 24.3050, "Ca" => 40.078, "Sr" => 87.62, "Ba" => 137.327}
  ATOM_VALENCES = {"H" => 1, "He" => 2, "B" => 3, "C" => 4, "N" => 5, "O" => 6, "F" => 7, "Ne" => 8, "Si" => 4, "P" => 5, "S" => 6, "Cl" => 7, "Ar" => 8, "As" => 5, "Se" => 6, "Br" => 7, "Kr" => 8, "Te" => 6, "I" => 7, "Xe" => 8, "At" => 7, "Rn" => 8, "Li" => 1, "Na" => 1, "K" => 1, "Rb" => 1, "Cs" => 1, "Be" => 2, "Mg" => 2, "Ca" => 2, "Sr" => 2, "Ba" => 2}
  PREFIXES = {"" => "mono", 2 => "di", 3 => "tri", 4 => "tetra", 5 => "penta", 6 => "hexa", 7 => "septa", 8 => "octa", 9 => "nona", 10 => "deca"}
  OXY_PREFIXES = {"" => "mon", 2 => "di", 3 => "tri", 4 => "tetr", 5 => "pent", 6 => "hex", 7 => "sept", 8 => "oct", 9 => "non", 10 => "dec"}
  IDE_ATOM_NAMES = {"H" => "hydride", "He" => "helium", "B" => "boride", "C" => "carbonide", "N" => "nitride", "O" => "oxide", "F" => "fluoride", "Ne" => "neon", "Si" => "silicide", "P" => "phosphide", "S" => "sulfide", "Cl" => "chloride", "Ar" => "argon", "As" => "arsenide", "Se" => "selenide", "Br" => "bromide", "Kr" => "krypton", "Te" => "telluride", "I" => "sodide", "Xe" => "xenon", "At" => "astatide", "Rn" => "radon"}
  BOND_TYPES = {1 => "single", 2 => "double", 3 => "triple", 4 => "quadruple"}

end

class String
  include Chemistry
  def contains_any_of?(array)
    array.each do |value|
      return true if self.include?(value)
    end
    false
  end

  # "N2O5" to [Atom.new("N"), Atom.new("N"), ] with 5 Atom.new("O"); Still cannot end with ")"
  def parse_as_symbol
    final = []
    characters = self.split("")
    i = 0

    while i < characters.length   #This whole long while loop is to take parentheses out; (CH2)2 -> CH2CH2, in all its varied forms;  it is still unrefined
      if characters[i] == ")" and characters[i+1] =~ /\d/
        holder = characters.join  #Next line alters characters - safety precaution
        block = characters.each_index{|x| characters[x] = characters[x].replace("") if x>i+1}.reject{|chr| chr == ""}  #Basically drop all after current index; aviods issues of array subtraction, likes array1 - ["3"] removes ALL instances of "3" in array1

        block = block.drop(1) until block.include?("(") == false  #Drop everything before "("; this leaves only the values between your parentheses
        block -= [")"] #Just a double check
        block.delete_at(block.length-1) if block.last =~ /\d/
        characters = holder.split("")
        str = ""
        (characters[i+1].to_i).times{|x| str += block.join}
        characters = characters.join.gsub!("(" + block.join + ")" + characters[i+1], str).split("")

        i = 0  #If something is changed, go back to the start of the string and recheck for parentheses
      end #End-if
      i += 1
    end #End-while
    
    characters = characters.join.gsub(/[(,)]/, "")  #Clean-up; Removes ( and ) in Zn(SCN) or CH(CCl4)

    characters.scan(/[A-Z][^A-Z]*/).each do |e| #The body of the method - letters as atom symbols, numbers as repeated atoms
      if e=~ /\d/
        e[/(\d)/, 1].to_i.times {|x| final << Atom.new(e.gsub(/\d/, ""))}

      else
        final << Atom.new(e)
      end #End-if
    end  #End-each

    final
  end

  # "Dinitrogen pentoxide" to same output as "N2O5"
  def get_symbols
    revised = self.downcase
    
    #The following is a systematic breakdown of covalent compound names into prefixes and atoms
    if revised.contains_any_of?(PREFIXES.values) #Replace each prefix with its number:  dinitrogen -> 2nitrogen
      #puts "contains a prefix"
      PREFIXES.values.each{|prefix| revised.gsub!(prefix, PREFIXES.invert[prefix].to_s) unless (revised.gsub(prefix, PREFIXES.invert[prefix].to_s).include?("xide") and revised.gsub(prefix, PREFIXES.invert[prefix].to_s).include?("oxide") == false)}  #The unless statement is to catch 'mono'xide
    end
    if revised.contains_any_of?(OXY_PREFIXES.values)  #Oxy prefixes; tetroxide -> 4oxide
      #puts "contains an oxy prefix"
      OXY_PREFIXES.values.each{|prefix| revised.gsub!(prefix, OXY_PREFIXES.invert[prefix].to_s)}
    end
    if revised.contains_any_of?(ATOM_NAMES.values)  #If there are any atoms names, replace each name with its symbol; 5oxide -> 5O
      #puts "contains an atom name"
      ATOM_NAMES.values.each{|value| revised.gsub!(value, ATOM_NAMES.invert[value])}
    end
    if revised.contains_any_of?(IDE_ATOM_NAMES.values) #Same, with ide names; 4phosphide -> 4P
      #puts "contains an ide atom name"
      IDE_ATOM_NAMES.values.each{|value| revised.gsub!(value, IDE_ATOM_NAMES.invert[value])}
    end

    final = ""
    #I stole this part - not really sure how it works
    revised.split(" ").each do |element|
      if element =~ /\d/
        element += element[/(\d)/,1]

        element = element[/.(.*)/m,1]
        final += element
      else
        final += element
      end #End-if

    end  #End-each
    final

  end
  
  #This is separated so that the information returned by .get_symbols can be found separately
  def parse_as_name  
    self.get_symbols.parse_as_symbol
  end

  #GREATEST METHOD EVER
  def is_symbolic_representation?
    self.match(/^([(]?[A-Z]*[a-z]?[0-9]?[)]?)*$/) != nil    
  end

end

#Goal: Atoms have names, symbols, weights, outer electron shell values (1..8), can bond based on octet rule; later, can ionize
class Atom
  include Chemistry
  attr_accessor :valenceshellnum, :bondcount, :charge
  attr_reader :weight, :name, :symbol
  #Make a new atom, make sure the element exists
  def initialize(symbol)
    begin
      raise ElementError, "What element is #{symbol}?" unless ATOM_NAMES.has_key?(symbol)
      @symbol = symbol
      @name = ATOM_NAMES[@symbol]
      @valenceshellnum = ATOM_VALENCES[@symbol]
      @weight = ATOM_WEIGHTS[@symbol]
      @bondcount = 0
      @charge = 0
    rescue ElementError
      puts $!
    end
  end

  #An atom is stable if it satisfies the octet rule; ionization will be a possibility later
  def stable?
    if self.symbol == "H" or self.symbol == "He"
      @valenceshellnum == 2 or @valenceshellnum == 0
    else
      @valenceshellnum == 8 or @valenceshellnum == 0
    end
  end

  def ionize!  #For atoms, not polyatomic ions; polyatomic ions ionize to fill octets not filled by bonding

  end

  #Bond an atom to another atom, increasing the valence number and bond count; passing bond_type allows for better error messages, used in other bond methods
  def s_bond!(other_atom, bond_type=1)
    begin

      if self.stable? == false and other_atom.stable? == false #neither atom has full octet
        self.valenceshellnum += 1
        other_atom.valenceshellnum += 1
        self.bondcount += 1
        other_atom.bondcount += 1
      elsif self.stable? and other_atom.stable? == false #self already has full octet
        raise ValenceShellError, "Cannot #{BOND_TYPES[bond_type]} bond: #{self.name} already has a full valence shell with #{self.bondcount} bonds"
      elsif self.stable? == false and other_atom.stable? #other_atom already has full octet
        raise ValenceShellError, "Cannot #{BOND_TYPES[bond_type]} bond: #{other_atom.name} already has a full valence shell with #{other_atom.bondcount} bonds"
      elsif self.stable? and other_atom.stable? #Both already have full octet
        raise ValenceShellError, "Cannot #{BOND_TYPES[bond_type]} bond: #{self.name} and #{other_atom.name} already have a full valence shells with #{self.bondcount} and #{other_atom.bondcount} bonds"
      end #End-if

    rescue ValenceShellError
      puts $!
    end
  end

  #Forms a double bond, which is the same as two single bonds
  def d_bond!(other_atom)
    2.times {|x| self.s_bond!(other_atom, x+1)}
  end

  #Triple bond
  def t_bond!(other_atom)
    3.times {|x| self.s_bond!(other_atom, x+1)}
  end

  #Quadruple bond
  def q_bond!(other_atom)
    4.times {|x| self.s_bond!(other_atom, x+1)}
  end

end

#Prereq - I suppose it must understand the Molecules are (trees, I'm going with trees) of Atoms
#Goal: Enter a compound name, get a chemical formula; enter a chemical formula, get a name as carbon dioxide <=> CO2; dinitrogen pentoxide <=> N2O5
class Molecule
  include Chemistry
  attr_accessor :symbolic, :name, :weight
  attr_reader :atoms
  def initialize(name)
    #Okay, here's a thought - If strict function formulae are require (Mean CH3COO not C2H3O2; CH3CH2Br not C2H5Br)
    @weight = 0
    if name.is_symbolic_representation?
      @atoms = name.parse_as_symbol
      @symbolic = name
      @atoms.each {|atom| @weight += atom.weight}
      @name = ""

    elsif name.is_symbolic_representation? == false
      @atoms = name.parse_as_name
      @name = name
      @atoms.each {|atom| @weight += atom.weight}
      @symbolic = @name.get_symbols

      
    end #End-if

  end

end

#Goal: Pretty much the same as Molecule, but as 2,3-dimethylbutane <=> CH3CH(CH3)CH(CH3)CH3; 3-methylpentane <=> CH3CH2CH(CH3)CH2CH3; 2-methylpentane <=> (CH3)2CH(CH2)2CH3
class OrganicMolecule < Molecule

end