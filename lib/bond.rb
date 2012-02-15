module Aims
  
=begin
A representation of a bond between two atoms
Bonds are not directional, so two bonds are equal if they bond the same atoms regradless of order
=end
class Bond

attr_accessor :atoms

	def initialize(atom1, atom2)
		self.atoms = [atom1, atom2]
	end

	def ==(bond)
		# Take the difference between the two sets
		diff1 = (self.atoms - bond.atoms)
		diff2 = (bond.atoms - self.atoms)
		# the lists are the same if both sets are empty
		diff1.empty? & diff2.empty?
	end

=begin
THe bond length is the distance between the two atoms
=end
	def length
		self.atoms[0].distance_to(self.atoms[1])
	end

=begin
Access the atoms in the bond
=end
	def [](i)
		atoms[i]
	end

end
end