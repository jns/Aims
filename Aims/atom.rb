module Aims
  
  class Atom
    @@lastid = 0

    attr_accessor :x, :y, :z, :id, :species, :constrain
    include Enumerable
    
    def initialize(x=nil, y=nil, z=nil, s=nil, c=Array.new)
      self.x = x
      self.y = y
      self.z = z
      self.species = s
      self.id = (@@lastid +=1)
      self.constrain = c
    end
=begin
  Two atoms are equal if their coordinates are equal and they are the same species
=end
	def ==(atom)
		(self.x == atom.x) & (self.y == atom.y) & (self.z == atom.z) & (self.species == atom.species)
	end

    def each
      [self.x, self.y, self.z].each{|i| yield i}
    end
    
    def [](i)
      case i
      when 0
        self.x
      when 1
        self.y
      when 2
        self.z
      end
    end
    
=begin
	Return the distance between to another atom
=end
	def distance_to(atom)
		Math.sqrt((self.x - atom.x)**2 + (self.y - atom.y)**2 + (self.z - atom.z)**2)
	end
	
    def copy
      Atom.new(self.x, self.y, self.z, self.species, self.constrain)
    end
    
    def displace(x,y,z)
      Atom.new(self.x+x, self.y+y, self.z+z, self.species, self.constrain)
    end

    def displace!(x,y,z)
      self.x += x
      self.y += y
      self.z += z
    end

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
    
    def rotate(mat)
      v = Vector[self.x, self.y, self.z]
      newv = mat*v
      Atom.new(newv[0], newv[1], newv[2], self.species, self.constrain)
    end
    
    def to_s
      "%s %16.6f %16.6f %16.6f" % [self.species, self.x, self.y, self.z]
    end
    
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