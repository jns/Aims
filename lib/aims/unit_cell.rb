module Aims

  class UnitCell

    include Vectorize
    include Enumerable
    
    attr_accessor :lattice_vectors
	  attr_accessor :bonds
	
    # clip_planes is an array of planes
    # The planes are defined with respect to the vectors millerX and millerZ
    # which define the miller indices corresponding to the x and z direction
    # Only atoms behind all the planes are visible
    # This can be used  to define a solid
    attr_accessor :clip_planes, :cart_to_miller, :miller_to_cart
    protected :clip_planes, :cart_to_miller, :miller_to_cart
    protected :clip_planes=, :cart_to_miller=, :miller_to_cart=


    def initialize(atoms, vectors = nil)  
      
      # Do some basic validation
      unless atoms.is_a? Array
        raise "Atoms must be an array!"
      end
      
      if atoms.empty? 
        raise "Atoms array is empty!"
      end
      
      atoms.each{|a| 
        unless a.is_a? Atom
          raise "Atoms array contains invalid object type #{a.class}!"
        end
      }
      # Ok. I'm satisfied
      self.atoms = atoms
      
      # Now check the lattice vectors
      if vectors
        self.lattice_vectors = vectors.collect{|v| 
          if v.is_a? Vector and v.size == 3
            v
          elsif v.is_a? Array and v.size == 3
            Vector.elements(v)
          else
            raise "Invalid lattice vector"
          end
        }
        unless self.lattice_vectors.size == 3
          raise "There must be 3 lattice vectors, not #{self.lattice_vectors.size}."
        end
      end
      
      @clip_planes = []      
      make_bonds
    end

	# Is the geometry empty
	def empty?
		@atoms.empty?
	end

    # Define the atoms in basis of this Unit Cell
    def atoms=(listOfAtoms)
      @atoms = listOfAtoms
      recache_visible_atoms
    end
    
    # Return the atoms in this unit cell.  Be default returns
    # only the visible atoms, but will return all of the atoms 
    # if called with (visibleOnly = false)
    def atoms(visibleOnly = true)
      if visibleOnly and 0 < @clip_planes.size
        @visibleAtoms
      else
        @atoms
      end
    end
    
	def make_bonds(bond_length = 4.0)
		# initialize an empty array
		self.bonds = Array.new
		
		# Make bonds between all atoms
		stack = atoms.dup
		
		atom1 = stack.pop
		while (not stack.empty?)
			stack.each{|atom2|
				b = Bond.new(atom1, atom2)
				self.bonds << b if b.length < bond_length
			}
			atom1 = stack.pop
		end
	end

  # Add a clip Plane to the unit cell
  def add_plane(aPlane, recache = true)
    self.clip_planes << aPlane
    recache_visible_atoms if recache
  end
	
    #  Add a clipping plane defined by the outward normal (h,k,l)
    # and a point on the plane (x,y,z)
    # Only points behind the plane will be kept.  Points in front
    # of the plane will be invisible.  They are not really gone, but 
    # can be returned by moving or removing the plane
    def add_plane_miller(h,k,l,x,y,z)
      normal = self.cartesian_from_miller(h, k, l)      
      self.add_plane_cartesian(normal[0], normal[1], normal[2], x, y, z)
    end
    
    #  Add a clipping plane defined by the outward normal (nx,ny,nz)
    # and a point on the plane (x,y,z)
    # Only points behind the plane will be kept.  Points in front
    # of the plane will be invisible.  They are not really gone, but 
    # can be returned by moving or removing the plane
    def add_plane_cartesian(nx, ny, nz, x, y, z)
      add_plane(Plane.new(nx, ny, nz, x, y, z), true)
    end
    
    # Clear all clip planes
    def clear_planes
      self.clip_planes.clear
      recache_visible_atoms
    end
    
    # remove a specific clip plane
    def remove_plane(aPlane)
      self.clip_planes.reject!{|p| p == aPlane}
      recache_visible_atoms
    end
    
    # Define the [h, k, l] vectors for each cartesian direction x, y, and z
    # There must be at least two vectors provided, and they must be orthogonal
    # The z vector, if provided, must also be orthogonal.  If it is not 
    # provided it will be calculated as the cross product of x and y
    # 
    # These vectors will populate a matrix A that will allow calculation of the
    # miller indices for any arbitrary vector in cartesian space, and 
    # the calculation of the cartesian space vector for any arbitrary miller index
    # @see miller_from_cartesion, cartesian_from_miller
    def set_miller_indices(x, y, z=nil)
      raise "Vectors must be orthogonal" unless 0 == dot(x,y)
      if z
        raise "Vectors must be orthogonal" unless 0 == dot(x,z)
      else
        z = cross(x, y)
      end
      self.cart_to_miller = Matrix[[x[0], y[0], z[0]], 
                               [x[1], y[1], z[1]],
                               [x[2], y[2], z[2]]]
      self.miller_to_cart = cart_to_miller.inverse
      return nil
    end

    # rotate the atoms about the z-axis so that the 
    # vector given by the miller index [h,k,l] is aligned
    # with the cartesian unit vector [1,0,0]
    # The vector [h,k,l] must be orthogonal to the defined
    # miller index of the z-direction
    def align_x(h,k,l)
      millerZ = self.cart_to_miller.column(2)
      unless 0 == dot([h,k,l], millerZ)
        raise "Specified vector [#{[h,k,l].join(',')}] is not orthogonal to z-axis [#{millerZ.to_a.join(',')}]"
      end

      # Define the current x axis and the new x-axis
      millerX = self.cart_to_miller.column(0)
      newX = Vector[h,k,l]

      # Find the angle between the current x direction and the new x-direction
      angle = acos(dot(newX*(1/newX.r), millerX*(1/millerX.r)))*180/Math::PI
      
      #Make the rotation in azimuth
      self.rotate(0, angle)

    end

    # Given a vector (x,y,z) in cartesian coordinates, 
    # return the miller-index corresponding to that vector's direction
    def miller_from_cartesian(x, y, z)
      if self.cart_to_miller
        self.cart_to_miller*Vector[x, y, z]
      else
        nil
      end
    end
    
    # Given a miller index (h,k,l), 
    # return a vector in cartesian coordinates pointing in that direction
    def cartesian_from_miller(h, k, l)

      if self.miller_to_cart
        self.miller_to_cart*Vector[h, k, l]
      else
        nil
      end
    end
    
    # Recompute the atoms that are behind all the clip planes
    def recache_visible_atoms

      plane_count = (@clip_planes ? @clip_planes.length : 0)
      return if plane_count == 0

      if @visibleAtoms
        @visibleAtoms.clear
      else
        @visibleAtoms = []
      end
      @atoms.each{|a|
        i = plane_count
        @clip_planes.each{|p|
          i = i-1 if 0 >= p.distance_to_point(a.x, a.y, a.z) 
        }
        @visibleAtoms << a if i == 0
      }
      
      make_bonds
    end
    
    def rotate(alt, az)
      altrad = Math::PI/180*alt
      azrad = Math::PI/180*az
      
      sinalt = Math::sin(altrad)
      cosalt = Math::cos(altrad)
      sinaz = Math::sin(azrad)
      cosaz = Math::cos(azrad)
      
      mat1 = Matrix[[cosaz, -sinaz, 0], [sinaz, cosaz, 0], [0, 0, 1]]
      mat2 = Matrix[[1.0, 0.0, 0.0], [0.0, cosalt, -sinalt], [0.0, sinalt, cosalt]]
      mat = mat1*mat2    
      newatoms = atoms.collect{|a|
        a.rotate(mat)
      }
      newvectors = lattice_vectors.collect{|v|
        mat*v
      }
      uc = UnitCell.new(newatoms, newvectors)
      uc.cart_to_miller = mat*self.cart_to_miller
      uc.miller_to_cart = uc.cart_to_miller.inverse

      return uc
    end
    
    def copy
      newAtoms = []
      newVecs = []
      @atoms.each{|a|
        newAtoms << a.copy
      }
      @lattice_vectors.each{|v|
        newVecs << v*1
      }
      uc = UnitCell.new(newAtoms, newVecs)
      uc.clip_planes = self.clip_planes
      uc.miller_to_cart = self.miller_to_cart
      uc.cart_to_miller = self.cart_to_miller
      uc.recache_visible_atoms
      return uc
    end

    def size
      atoms.size
    end
    
    def bounding_box(visible_only = true)
      maxX = atoms(visible_only).first.x
      maxY = atoms(visible_only).first.y
      maxZ = atoms(visible_only).first.z
      minX = maxX
      minY = maxY
      minZ = maxZ
      
      atoms(visible_only).each{|a|
        if a.x > maxX
          maxX = a.x
        elsif a.x < minX
          minX = a.x
        end
        
        if a.y > maxY
          maxY = a.y
        elsif a.y < minY
          minY = a.y
        end
        
        if a.z > maxZ
          maxZ = a.z
        elsif a.z < minZ
          minZ = a.z
        end
      }

      [Atom.new(maxX, maxY, maxZ), Atom.new(minX, minY, minZ)]
    end
    
    def center(visible_only = true)
      bounds = bounding_box(visible_only)
      x = (bounds[0].x + bounds[1].x)/2.0
      y = (bounds[0].y + bounds[1].y)/2.0
      z = (bounds[0].z + bounds[1].z)/2.0
      return Atom.new(x,y,z)
    end
    
    def each
      self.atoms.each{|a|
        yield a
      }
    end
    
    def remove_atom(atom)
      atoms.reject!{|a|
       a.id == atom.id 
      }
      # Force a rehash of nearest-neighbor tree
      @tree = nil
    end
    
    def [](index)
      atoms[index]
    end
    
=begin
Return a new unit cell with all the atoms displaced by the amount x,y,z
=end
    def displace(x,y,z)
      UnitCell.new(atoms.collect{|a|
        a.displace(x,y,z)
      }, self.lattice_vectors)
      #TODO copy miller indices
    end

=begin
Repeat a unit cell nx,ny,nz times in the directions 
of the lattice vectors.
Negative values of nx,ny or nz results in displacement in the
negative direction of the lattice vectors
=end
    def repeat(nx=1, ny=1, nz=1)
      
      raise "Not a periodic system." if self.lattice_vectors.nil?
      
      u = self.copy
      v1 = self.lattice_vectors[0]
      v2 = self.lattice_vectors[1]
      v3 = self.lattice_vectors[2]

      nx_sign = (0 < nx) ? 1 : -1
      ny_sign = (0 < ny) ? 1 : -1
      nz_sign = (0 < nz) ? 1 : -1
      
      new_atoms = []
      nx.abs.times do |i|        
        ny.abs.times do |j|
          nz.abs.times do |k|
            new_atoms << self.displace(nx_sign*i*v1[0] + ny_sign*j*v2[0] + nz_sign*k*v3[0], 
                                       nx_sign*i*v1[1] + ny_sign*j*v2[1] + nz_sign*k*v3[1], 
                                       nx_sign*i*v1[2] + ny_sign*j*v2[2] + nz_sign*k*v3[2]).atoms
          end
        end
      end

      u.atoms = new_atoms.flatten
      u.lattice_vectors = [Vector[nx.abs*v1[0], nx.abs*v1[1], nx.abs*v1[2]], 
                           Vector[ny.abs*v2[0], ny.abs*v2[1], ny.abs*v2[2]], 
                           Vector[nz.abs*v3[0], nz.abs*v3[1], nz.abs*v3[2]]]
      u.make_bonds 
      return u     
    end

=begin
Concatenate the atoms from another unit cell to this unit cell.
Currently does no validation on lattice vectors on miller vectors
=end    
    def <<(aUnitCell)
      self.atoms.concat(aUnitCell.atoms)
	  self.make_bonds
	  return self
    end
    
    def to_s
      self.atoms.collect{|a| a.to_s}.join("\n")
    end
    
    def format_geometry_in
      output = ""
      if self.lattice_vectors
        output << self.lattice_vectors.collect{|v| "lattice_vector #{v[0]} #{v[1]} #{v[2]}"}.join("\n")
        output << "\n"
      end
      output << self.atoms.collect{|a| a.format_geometry_in}.join("\n")

      output
    end

=begin
return a string in xyz format
=end
    def format_xyz
      output = self.atoms.size.to_s + "\n"
      output << "Aims Geometry \n"
      self.atoms.each{ |a| 
        output << [a.species, a.x.to_s, a.y.to_s, a.z.to_s].join("\t") + "\n"
      }
      output
    end
    
    # Find the difference between this cell and another cell
    # Return a cell with Pseudo-Atoms whose positions are really the differences
    def delta(aCell)
      raise "Cells do not have the same number of atoms" unless self.atoms.size == aCell.atoms.size

      pseudo_atoms = []
      self.atoms.size.times {|i|
        a1 = self.atoms[i]
        a2 = aCell.atoms[i]
        raise "Species do not match" unless a1.species == a2.species
        a = Atom.new
        a.species = a1.species
        a.x = a1.x - a2.x
        a.y = a1.y - a2.y
        a.z = a1.z - a2.z
        pseudo_atoms << a
      }
      UnitCell.new(pseudo_atoms)
    end
    

=begin
  Move all atoms inside the primitive volume defined by the
  six planes of the lattice vectors
=end    
    def correct
      
      # Hash for storing bounding planes and the out-of-plane vector
      # by which each atom will be displaced to move it into the primitive volume
      # key = bounding plane
      # value = out-of-plane lattice vector used to displace atoms
      planes_vecs = {}

      # Define the primitive volume as six planes by 
      # finding the normal to each pair of lattice vectors
      # and making a plane with this normal that includes
      #   1. The point at the head of the third lattice vector 
      #      and pointing in the positive direction
      #
      #   2. The point at the tail of the third lattice vector
      #      and pointing in the negative direction  
      #
      (0..2).each do |i|
        out_vector = lattice_vectors[i] # The out of plane vector
        plane_vectors = lattice_vectors.reject{|v| v == out_vector} # The in plane vectors
        norm = cross(plane_vectors[0], plane_vectors[1])
        # if the norm has a component in the direction of the out of plane vector
        # then use the head of the out-of-plane vector as the intersection point
        # otherwise use the tail (the origin)
        if 0 < dot(norm, out_vector)
          # First plane is in direction of norm and intersects head
          # Displace vector is +1 if plane intersects tail and -1 if plane intersects head.
          planes_vecs[Plane.new(norm[0], norm[1], norm[2], out_vector[0], out_vector[1], out_vector[2])] = out_vector*(-1)
          # Second plane is opposite direction of norm and intersects tail
          planes_vecs[Plane.new(-norm[0], -norm[1], -norm[2], 0, 0, 0)] = out_vector*1
        else
          # First plane is in opposite direction of norm and intersects head
          planes_vecs[Plane.new(-norm[0], -norm[1], -norm[2], out_vector[0], out_vector[1], out_vector[2])] = out_vector*(-1)
          # Second plane is in direction of norm and intersects tail
          planes_vecs[Plane.new(norm[0], norm[1], norm[2], 0, 0, 0)] = out_vector*1
        end
      end

      # Make a coyp of the unit cell
      new_unit_cell = self.copy

      # Move each atom behind all the planes
      new_unit_cell.atoms(false).each do |atom|
        planes_vecs.each_pair do |p, v|
          if p.distance_to_point(0,0,0) == 0
            # If the plane intersects the origin then 
            # move atoms not on the plane (inequality)
            while p.distance_to_point(atom.x, atom.y, atom.z) > 0
              atom.displace!(v[0], v[1], v[2])
            end
          else
            # Move atoms that lie on the plane if the plane doesn't intersect the origin
            while p.distance_to_point(atom.x, atom.y, atom.z) >= 0
              atom.displace!(v[0], v[1], v[2])
            end
          end
        end
      end
      new_unit_cell.atoms.uniq!
      new_unit_cell.make_bonds
      return new_unit_cell

    end

  end
end
