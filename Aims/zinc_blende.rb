require 'mathn'

module Aims

# Factory class for generating slabs for various crystal surfaces
# of specified thickness and with the specified vacuum. 
  class ZincBlende

    include Math
    
    attr_accessor :lattice_const, :cation, :anion

    # Initialize the zinc-blende UnitCell
    # cation and anion are the  atomic 
    # species occupying the two different sub-lattices.
    # lattice_const specifies the lattice constant
    def initialize(cation, anion, lattice_const)
      self.lattice_const = lattice_const
      self.cation = cation
      self.anion = anion
    end
    
    def get_bulk
      b = 0.25*self.lattice_const
      a1 = Atom.new(0, 0, 0, self.cation)
      a2 = Atom.new(b, b, b, self.anion)
      
      v1 = Vector[0.5, 0.5, 0.0]*self.lattice_const
      v2 = Vector[0.5, 0.0, 0.5]*self.lattice_const
      v3 = Vector[0.0, 0.5, 0.5]*self.lattice_const
      zb = UnitCell.new([a1, a2], [v1, v2, v3])
      
      millerx = [1, 0, 0]
      millery = [0, 1, 0]
      millerz = [0, 0, 1]
      
      zb.set_miller_indices(millerx, millery, millerz)
      
      return zb
    end

    # Return a unit cell for a slab of 001 
    # Specify the number of atomic monolayers 
    # and the vacuum thickness in angstrom
    def get_001_surface(monolayers, vacuum)
      
    end

    # Return a unit cell for a slab of 111A (anion terminated)
    # specify the number of atomic monolayers 
    # and the vacuum thickness in angstrom
    def get_111A_surface(monolayers, vacuum)

    end

    # Return a unit cell for a slab of 111B (cation terminated)
    # specify the number of atomic monolayers 
    # and the vacuum thickness in angstrom
    def get_111B_surface(monolayers, vacuum)

      # The atoms on a FCC 
      as1 = Atom.new(0.0, 0.0, 0.0, self.cation)
      ga1 = Atom.new(0.0, 0.0, -sqrt(3)/4*self.lattice_const, self.anion)

      # The lattice Vectors
      v1 = Vector[0.5*sqrt(2), 0.0, 0.0]*self.lattice_const
      v2 = Vector[sqrt(2)*0.25, sqrt(6)*0.25, 0.0]*self.lattice_const
      v3 = Vector[sqrt(2)*0.25, sqrt(2.0/3.0)*0.25, -1*sqrt(4.0/3.0)*0.5]*self.lattice_const

      # The unit cell
      zb = UnitCell.new([as1, ga1], [v1, v2, v3])

      # The Miller Indices
      millerX = [-1, 1, 0]  # Orientation of the crystal pointing in the cartesian +x axis
      millerY = [1, 1, -2]  # Orientation of the crystal pointing in the cartesian +y axis
      millerZ = [-1, -1, -1] # Orientation of the crystal pointing in the cartesian +z axis

      zb.set_miller_indices(millerX, millerY, millerZ)

      # Repeat the unit cell and add vacuum
      if 0 < vacuum 
        layers = 9
        zb = zb.repeat(1,1,monolayers)
        bilayerSep = v3[2]

        zb.lattice_vectors[2] = Vector[0, 0, monolayers*bilayerSep + vacuum]

        # Constrain the bottom 2 layers
        # zb.atoms.each{|a|
        #   if (a.z < bilayerSep*2.0)
        #     a.constrain = ".true."
        #   end
        # }
      end
      
      zb
    end

    # return a unit cell for a slab of 112 
    # specify the number of atomic monolayers and the vacuum thickness in angstrom
    def get_112_surface(monolayers, vacuum=0)
      atom1 = Atom.new(0,0,0,self.cation)
      atom2 = Atom.new(self.lattice_const*sqrt(3)/2, 0, 0, self.anion)
      
      v1 = Vector[sqrt(3), 0, 0]*self.lattice_const
      v2 = Vector[0, sqrt(2)/2, 0]*self.lattice_const
      v3 = Vector[1/sqrt(3), 1/(sqrt(3)*2), -1/(sqrt(3)*2)]*self.lattice_const
      
      millerX = Vector[1, 1, -2];
      millerY = Vector[-1, 1, 0];
      millerZ = Vector[-1, -1, -1]

      # The unit cell
       zb = UnitCell.new([atom1, atom2], [v1, v2, v3])
       zb.set_miller_indices(millerX, millerY, millerZ)

       # Repeat the unit cell
       zb = zb.repeat(1,1,monolayers)


       if 0 < vacuum
         # Add vacuum
         monolayerSep = v3[2]
         zb.lattice_vectors[2] = Vector[0, 0, (monolayers*monolayerSep).abs + vacuum.to_f]
         # Move everything into a nice tidy unit cell. 
         zb = zb.correct
       end

       # # Constrain the bottom 2 layers
       # zb.atoms.each{|a|
       #   if (a.z < monolayerSep*2)
       #     a.constrain = ".true."
       #   end
       # }


       # Return the completed unit cell
       return zb
    end
    
    # Return a unit cell for a slab of 110
    # specify the number of atomic monolayers 
    # and the vacuum thickness in angstrom
    def get_110_surface(monolayers, vacuum=0)

      # The atoms on a FCC 
      atom1 = Atom.new(0,0,0,self.cation)
      atom2 = Atom.new(self.lattice_const*1/(2*sqrt(2)), self.lattice_const*0.25, 0.0, self.anion)

      # The lattice Vectors
      v1 = Vector[1/sqrt(2), 0.0, 0.0]*self.lattice_const
      v2 = Vector[0.0, 1.0, 0.0]*self.lattice_const
      v3 = Vector[1/(2*sqrt(2)), -0.5, 1/(2*sqrt(2))]*self.lattice_const

      # The miller indices for each primitive cartesian direction
      millerX = Vector[1, -1, 0]
      millerY = Vector[0, 0, 1]
      millerZ = Vector[1, 1, 0]

      # The unit cell
      zb = UnitCell.new([atom1, atom2], [v1, v2, v3])
      zb.set_miller_indices(millerX, millerY, millerZ)

      # Repeat the unit cell
      zb = zb.repeat(1,1,monolayers)


      if 0 < vacuum
        # Add vacuum
        monolayerSep = v3[2]
        zb.lattice_vectors[2] = Vector[0, 0, monolayers*monolayerSep + vacuum.to_f]
        # Move everything into a nice tidy unit cell. 
        zb = zb.correct
      end

      # # Constrain the bottom 2 layers
      # zb.atoms.each{|a|
      #   if (a.z < monolayerSep*2)
      #     a.constrain = ".true."
      #   end
      # }

      
      # Return the completed unit cell
      return zb
    end
  end

end
