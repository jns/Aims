module Aims
  
  class Atom
    # The last id assigned to an atom
    @@lastid = 0

    # The x coordinate of the atom in angstrom
    attr_accessor :x
    # The y coordinate of the atom in angstrom
    attr_accessor :y
    # The z coordinate of the atom in angstrom
    attr_accessor :z
    # The +id+ of this atom. Every atom has a unique id
    attr_accessor :id
    # The species of this atom
    attr_accessor :species
    # The relaxation constraints of this atom
    attr_accessor :constrain
    # Two atoms are equal if their coordinates are the same to this precision
    attr_accessor :precision
    
    include Enumerable
    
    # Create an atom of the specified species at the given coordinates
    # * +x+ The x coordinate of the atom in angstrom
    # * +y+ The y coordinate of the atom in angstrom
    # * +z+ The z coordinate of the atom in angstrom
    # * +s+ The atomic species ex. "C", "Si", "S", etc. (can be nil)
    # * +c+ The relaxation constraints. valid values are TRUE, FALSE, ".true.", ".false.", "x", "y", "z" or %w(x y z)
    def initialize(x=nil, y=nil, z=nil, s=nil, c=Array.new)
      self.x = x
      self.y = y
      self.z = z
      self.species = s
      self.precision = 0.0001
      self.id = (@@lastid +=1)
      self.constrain = c
    end
    
    # A boolean value, 
    # True if the atom has relaxation constraints
    def constrained?
      if self.constrain
        if self.constrain == true
          true
        elsif self.constrain.is_a? String
          true
        elsif self.constrain.is_a? Array and not self.constrain.empty?
          true
        else
          false
        end
      else
        false
      end
    end
    
  # Two atoms are equal if their coordinates are equal and they are the same species
	def ==(atom)
		((self.x-atom.x).abs < self.precision) & 
		((self.y-atom.y).abs < self.precision) & 
		((self.z-atom.z).abs < self.precision) & 
		(self.species == atom.species)
	end
  alias_method :eql?, :==
    
  # Implementation for Hash equality testing
  def hash
    (self.x*self.y + self.z).abs.ceil
  end
  
  # Enumerate over each coordinate (x,y,z)
  def each
    [self.x, self.y, self.z].each{|i| yield i}
  end
    
  # Index into the Atom's coordinates (x,y,z)
  def [](i)
    case i
    when 0
      self.x
    when 1
      self.y
    when 2
      self.z
    else
      raise "Index Out of Bounds"
    end
  end
    
  # Return the distance to another atom
	def distance_to(atom)
		Math.sqrt((self.x - atom.x)**2 + (self.y - atom.y)**2 + (self.z - atom.z)**2)
	end
	
	# A deep copy of the atom
    def copy
      Atom.new(self.x, self.y, self.z, self.species, self.constrain)
    end
    
    # Return a new atom with the same species and relaxation constraints
    # but with coordinates displaced by +x+, +y+, +z+
    def displace(x,y,z)
      Atom.new(self.x+x, self.y+y, self.z+z, self.species, self.constrain)
    end

    # Displace this atom in place
    def displace!(x,y,z)
      self.x += x
      self.y += y
      self.z += z
    end

    # Return an atom rotated about the z-axis using the origin as the center-point.
    # * +angle+ Is the amount to rotate in degrees (or it can respond to :sin and :cos) 
    def rotate_Z(angle)
      sinA = if angle.respond_to? :sin
        angle.sine
      else
        Math.sin(angle*Math::PI/180)
      end
      cosA = if angle.respond_to? :cos
        angle.cos
      else
        Math.cos(angle*Math::PI/180)
      end
      
      mat = Matrix[[cosA, -1*sinA, 0],[sinA, cosA, 0], [0,0,1]]
      rotate(mat) 
    end
    
    # Return an atom rotated about the x-axis using the origin as the center-point.
    # * +angle+ Is the amount to rotate in degrees (or it can respond to :sin and :cos) 
    def rotate_X(angle)
      sinA = if angle.respond_to? :sin
        angle.sine
      else
        Math.sin(angle*Math::PI/180)
      end
      cosA = if angle.respond_to? :cos
        angle.cos
      else
        Math.cos(angle*Math::PI/180)
      end
      mat = Matrix[[1, 0, 0], [0, cosA, -1*sinA],[0, sinA, cosA]]
      rotate(mat)       
    end
    
    # Return an atom rotated about the y-axis using the origin as the center-point.
    # * +angle+ Is the amount to rotate in degrees (or it can respond to :sin and :cos) 
    def rotate_Y(angle)
      sinA = if angle.respond_to? :sin
        angle.sine
      else
        Math.sin(angle*Math::PI/180)
      end
      cosA = if angle.respond_to? :cos
        angle.cos
      else
        Math.cos(angle*Math::PI/180)
      end
      mat = Matrix[[cosA, 0, -1*sinA],[0, 1, 0], [sinA, 0, cosA]]
      rotate(mat)             
    end

    # Return a new rotated atom about the origin using the given 3x3 Math::Matrix.
    def rotate(mat)
      v = Vector[self.x, self.y, self.z]
      newv = mat*v
      Atom.new(newv[0], newv[1], newv[2], self.species, self.constrain)
    end
    
    # Print a string representation of this atom
    def to_s
      "%s %16.6f %16.6f %16.6f" % [self.species, self.x, self.y, self.z]
    end
    
    # Print a string representation of this atom formatted in the 
    # geometry.in format used by Aims
    def format_geometry_in
      line = "atom %16.6f %16.6f %16.6f %s" % [self.x, self.y, self.z, self.species]
      if self.constrain
        if self.constrain == true
          line << "\nconstrain_relaxation .true."
        elsif self.constrain.is_a? String
          line << "\nconstrain_relaxation #{self.constrain}"
        elsif self.constrain.is_a? Array and not self.constrain.empty?
          self.constrain.each{|c|
            line << "\nconstrain_relaxation #{c}"            
          }
          line << "\n"
        end
      end
      line
    end
  end
end