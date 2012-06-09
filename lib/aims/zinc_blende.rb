require 'mathn'

module Aims

# Factory class for generating slabs for various crystal surfaces
# of specified thickness and with the specified vacuum. 
# Example:
#   zb = ZincBlende.new("Ga", "As", 5.65)
#   zb.get_bulk -> # returns the a unit cell with two atoms
#   zb.get_111B_surface(7, 20, 3) # returns a slab of 111B 7 monolayers thick, with 20 angstrom of vacuum, and the bottom 3 layers constrained
  class ZincBlende

    include Vectorize
    include Math
    
    attr_accessor :lattice_const, :cation, :anion

    # Initialize the zinc-blende Geometry
    # cation and anion are the  atomic 
    # species occupying the two different sub-lattices.
    # lattice_const specifies the lattice constant
    def initialize(cation, anion, lattice_const)
      self.lattice_const = lattice_const
      self.cation = cation
      self.anion = anion
    end
    
    # Return the traditional unit cell of bulk zinc blende
    def get_bulk
      b = 0.25*self.lattice_const
      a1 = Atom.new(0, 0, 0, self.cation)
      a2 = Atom.new(b, b, b, self.anion)
      
      v1 = Vector[0.5, 0.5, 0.0]*self.lattice_const
      v2 = Vector[0.5, 0.0, 0.5]*self.lattice_const
      v3 = Vector[0.0, 0.5, 0.5]*self.lattice_const
      zb = Geometry.new([a1, a2], [v1, v2, v3])
      
      millerx = [1, 0, 0]
      millery = [0, 1, 0]
      millerz = [0, 0, 1]
      
      zb.set_miller_indices(millerx, millery, millerz)
      
      return zb
    end

    # Fill the given volume with atoms
    def fill_volume(volume)
      
      # First fill a cube that bounds the volume
      max = volume.max_point
      min = volume.min_point
      
      dx = max[0] - min[0]
      dy = max[1] - min[1]
      dz = max[2] - min[2]
      
      bulk = get_bulk
      
      # This inverse matrix gives the number of repetitions 
      m = Matrix[[dx,0,0], [0,dy,0], [0,0,dz]]
      v = Matrix[bulk.lattice_vectors[0].to_a, 
                 bulk.lattice_vectors[1].to_a,
                 bulk.lattice_vectors[2].to_a]
      rep_mat = m*(v.inverse)
      
      # The only way I can figure out how to do this for an 
      # arbitrary set of lattice vectors is to fill the volume
      # out along each edge of the super-cube and then eliminate duplicates
      atoms = []
      3.times do |i|
        # this vector is the number of repetitions in the unit cell
        # to fill the volume out along the i-th edge of the super-cube
        n_repeat = rep_mat.row(i)
        
        # Give the proper sign to the repeat
        nx = (n_repeat[0] < 0) ? n_repeat[0].floor-1 : n_repeat[0].ceil+1
        ny = (n_repeat[1] < 0) ? n_repeat[1].floor-1 : n_repeat[1].ceil+1
        nz = (n_repeat[2] < 0) ? n_repeat[2].floor-1 : n_repeat[2].ceil+1
        
        atoms += bulk.repeat(nx, ny, nz).atoms.find_all{|a| volume.contains_point(a.x, a.y, a.z)}
      end
      Geometry.new(atoms.uniq)
    end

    # Return a unit cell for a slab of 001 
    # Specify the number of atomic monolayers,
    # the vacuum thickness in angstrom,
    # and the number of layers to constrain at the base of the slab
    def get_001_surface(monolayers, vacuum, constrain_layers = 0)
      anion = Atom.new(0,0,0,self.cation)
      cation = Atom.new(0.25*self.lattice_const, 0.25*self.lattice_const, 0.25*self.lattice_const, self.anion)
      v1 = Vector[0.5, 0.5, 0]*self.lattice_const
      v2 = Vector[-0.5,0.5,0]*self.lattice_const
      v3 = Vector[0.5, 0, 0.5]*self.lattice_const
      
      zb = Geometry.new([anion, cation], [v1,v2,v3])
      millerX = [1,0,0]
      millerY = [0,1,0]
      millerZ = [0,0,1]
      zb.set_miller_indices(millerX, millerY, millerZ)
      
      # Repeat the unit cell.  The unit cell is a bi-layer so divide by 2
      zb = zb.repeat(1,1,(monolayers/2).ceil)

      if 0 < vacuum
        # Add vacuum
        monolayerSep = v3[2]/2
        zb.lattice_vectors[2] = Vector[0, 0, (monolayers-1)*monolayerSep.abs + vacuum.to_f]
        # Move everything into a nice tidy unit cell. 
        zb = zb.correct
      end

      minZ = zb.atoms.min{|a,b| a.z <=> b.z}.z
      
      # Reject the top layer of atoms if an odd number of monolayers was requested.
      # This is necessary because the primitive cell is a bilayer
      zb.atoms.reject! {|a| 
        a.z >= (minZ + monolayerSep.abs*monolayers)
      }
      
      # Constrain the bottom layers
      zb.atoms.each{|a|
        if (a.z < minZ + monolayerSep.abs*constrain_layers)
          a.constrain = ".true."
        end
      }
      
      # Return the completed unit cell
      return zb
    end
    alias_method :get_100_surface, :get_001_surface
    alias_method :get_010_surface, :get_001_surface

    # Return a unit cell for a slab of 111A (anion terminated)
    # specify the number of atomic monolayers,
    # the vacuum thickness in angstrom,
    # and the number of layers to constrain at the base of the slab
    def get_111A_surface(monolayers, vacuum, constrain_layers = 0)
      # get the 111B surface, then reflect it about z=0
      get_111_surface("A", monolayers, vacuum, constrain_layers)
    end

    # Return a unit cell for a slab of 111A (cation terminated)
    # specify the number of atomic monolayers,
    # the vacuum thickness in angstrom,
    # and the number of layers to constrain at the base of the slab
    def get_111B_surface(monolayers, vacuum, constrain_layers = 0)
      get_111_surface("B", monolayers, vacuum, constrain_layers)
    end

    # Return a unit cell for a slab of 111
    # dir is either "A" or "B" for the cation or anion terminated slab
    # specify the number of atomic monolayers 
    # and the vacuum thickness in angstrom
    def get_111_surface(dir, monolayers, vacuum, constrain_layers = 0)

      if dir == "A"
        top_atom = self.anion
        bot_atom = self.cation
      elsif dir == "B"
        top_atom = self.cation
        bot_atom = self.anion
      else
        raise "Direction must be either A or B"
      end
      
      # The atoms on a FCC 
      as1 = Atom.new(0.0, 0.0, 0.0, top_atom)
      ga1 = Atom.new(0.0, 0.0, -sqrt(3)/4*self.lattice_const, bot_atom)

      # The lattice Vectors
      v1 = Vector[0.5*sqrt(2), 0.0, 0.0]*self.lattice_const
      v2 = Vector[sqrt(2)*0.25, sqrt(6)*0.25, 0.0]*self.lattice_const
      v3 = Vector[sqrt(2)*0.25, sqrt(2.0/3.0)*0.25, -1*sqrt(4.0/3.0)*0.5]*self.lattice_const

      # The unit cell
      zb = Geometry.new([as1, ga1], [v1, v2, v3])

      # The Miller Indices
      millerX = [-1, 1, 0]  # Orientation of the crystal pointing in the cartesian +x axis
      millerY = [1, 1, -2]  # Orientation of the crystal pointing in the cartesian +y axis
      millerZ = [-1, -1, -1] # Orientation of the crystal pointing in the cartesian +z axis

      zb.set_miller_indices(millerX, millerY, millerZ)

      # Repeat the unit cell and add vacuum
      if 0 < vacuum 
        # We actually repeat the unit cell monolayers+1 times because
        # I will strip off the top and bottom atoms to make the proper surface
        zb = zb.repeat(1,1,monolayers+1)
        
        bilayerSep = v3[2]
        zb.lattice_vectors[2] = Vector[0, 0, (monolayers-1)*(bilayerSep.abs) + vacuum]

        # Strip off the top and bottom atom
        minZ = zb.atoms.min{|a,b| a.z <=> b.z}.z
        maxZ = zb.atoms.max{|a,b| a.z <=> b.z}.z

        zb.atoms.reject!{|a| a.z == maxZ}
        zb.atoms.reject!{|a| a.z == minZ}

        # Constrain the bottom layers if requested
        if 0 < constrain_layers
          # get the min again because we removed the atoms at minZ above
          minZ = zb.atoms.min{|a,b| a.z <=> b.z}.z
          constrain_below = minZ + bilayerSep.abs*constrain_layers
          zb.atoms.each{|a|
            if (a.z < constrain_below)
              a.constrain = ".true."
            end
          }
        end
      end
      
      zb
    end

    # return a unit cell for a slab of 112 
    # specify the number of atomic monolayers and the vacuum thickness in angstrom
    def get_112_surface(monolayers, vacuum=0, constrain_layers = 0)
      atom1 = Atom.new(0,0,0,self.cation)
      atom2 = Atom.new(self.lattice_const*sqrt(3)/2, 0, 0, self.anion)
      
      v1 = Vector[sqrt(3), 0, 0]*self.lattice_const
      v2 = Vector[0, sqrt(2)/2, 0]*self.lattice_const
      v3 = Vector[1/sqrt(3), 1/(sqrt(3)*2), -1/(sqrt(3)*2)]*self.lattice_const
      
      millerX = Vector[1, 1, -2];
      millerY = Vector[-1, 1, 0];
      millerZ = Vector[-1, -1, -1]

      # The unit cell
       zb = Geometry.new([atom1, atom2], [v1, v2, v3])
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
    def get_110_surface(monolayers, vacuum=0, constrain_layers = 0)

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
      zb = Geometry.new([atom1, atom2], [v1, v2, v3])
      zb.set_miller_indices(millerX, millerY, millerZ)

      # Repeat the unit cell
      zb = zb.repeat(1,1,monolayers)


      monolayerSep = v3[2]
      if 0 < vacuum
        # Add vacuum
        zb.lattice_vectors[2] = Vector[0, 0, (monolayers-1)*monolayerSep.abs + vacuum.to_f]
        # Move everything into a nice tidy unit cell. 
        zb = zb.correct
      end

      # # Constrain the bottom layers
      zb.atoms.each{|a|
        if (a.z < monolayerSep*constrain_layers)
          a.constrain = ".true."
        end
      }

      
      # Return the completed unit cell
      return zb
    end
  end

end
