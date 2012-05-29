
require 'mathn'

module Aims
  class Wurtzite
    
    include Math
    
    attr_accessor :lattice_const, :cation, :anion

    # Initialize the wurtzite Geometry
    # cation and anion are the  atomic 
    # species occupying the two different sub-lattices.
    # lattice_const specifies the lattice constant
    def initialize(cation, anion, lattice_const)
      self.lattice_const = lattice_const
      self.cation = cation
      self.anion = anion
    end
    
    def get_bulk
      # The lattice constant
      a = (ARGV[0] || 4.0).to_f
      c = 1.63299*a # sqrt(8/3)a
      #a = 6.19231 # Strained GaSb, Unstrained GaSb is 6.22
      #a = 5.75 # Unstrained GaAs

      # The atoms on a HCP
      as1 = Atom.new(0,0,0,'As')
      as2 = as1.displace(0.5*a, 0.33333*a, 0.5*c)

      ga1 = Atom.new(0.0, 0.0, c*0.386, 'Ga')
      ga2 = ga1.displace(0.5*a,0.33333*a, 0.5*c)

      # The lattice Vectors
      v1 = Vector[1.0, 0.0, 0.0]*a
      v2 = Vector[0.5, 0.5*sqrt(3), 0.0]*a
      v3 = Vector[0.0, 0.0, 1.0]*c

      # The unit cell
      wz = Geometry.new([as1,ga1,as2,ga2], [v1, v2, v3])
      
      # wz.set_miller_indices(millerX, millerY, millerZ)
      return wz
    end
  end
end
   