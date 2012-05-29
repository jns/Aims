module Aims
  
# A representation of a bond between two atoms
# Bonds are not directional, so two bonds are equal if they bond the same atoms regradless of order
class Bond

attr_accessor :atoms

  # Initialize a bond between two atoms
	def initialize(atom1, atom2)
		self.atoms = [atom1, atom2]
	end

  # Two bonds are equal iff the set of atoms
  # in both bonds is the same.
	def ==(bond)
		# Take the difference between the two sets
		diff1 = (self.atoms - bond.atoms)
		diff2 = (bond.atoms - self.atoms)
		# the lists are the same if both sets are empty
		diff1.empty? & diff2.empty?
	end

  # Two bonds are eql? iff the set of atoms
  # in both bonds is the same.
  def eql?(bond)
    self == bond
  end

  # Implementation of hash for equality testing
  def hash
    self.atoms.hash
  end

  # The bond length is the distance between the two atoms
	def length
		self.atoms[0].distance_to(self.atoms[1])
	end

  # Access the atoms in the bond
	def [](i)
		atoms[i]
	end

end
end