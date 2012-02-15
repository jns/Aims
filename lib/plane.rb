module Aims
  
  class Plane
=begin
 A Class representing a plane.
 Internally stored in Hessian Normal Form.
 0 = Ax + By + Cz - D
 [A,B,C] is the plane normal.
 D is the distance from the origin.
=end
    attr_reader :a, :b, :c, :d
  
    # Initialize this plane with the normal and a
    # point (x,y,z) on the plane
    def initialize(a, b, c, x, y,z)
      @a = a
      @b = b
      @c = c
      if (@a == 0 and @b == 0 and @c == 0)
        raise "Invalid definition of plane."
      end
      @d = a*x + b*y + c*z
    end

    # return some arbitrary point on the plane
    def any_point_on_plane
      
      unless (@c == 0)
        return Vector[0, 0, @d/@c]
      end
      
      unless (@b == 0)
        return Vector[0, @d/@b, 0]
      end
      
      unless (@a == 0)
        return Vector[@d/@a, 0, 0]
      end
      
      # Actually if we get to this point, the plane undetermined and all of R3 satisfies the definition
      return Vector[0,0,0]
    end
    
    # The distance to a point (x,y,z) is
    # D - Ax - By - Cz
    def distance_to_point(x, y, z)
      a*x + b*y + c*z - d
    end
    
    # Return the unit normal
    def unit_normal
      v = Vector[@a, @b, @c]
      v*(1/v.r)
    end
    
    # The equation for the interstion of a ray and a plane
    #  
    def intersection_with_ray(a, b)
      
    end
    
    # Displace this plane in the direction of its normal
    def displace_along_normal(distance)
      @d += distance
    end
    
    # Two planes are equal if there ABCD parametersequal
    def ==(aPlane)
      @a == aPlane.a and @b == aPlane.b and @c == aPlane.c and @d == aPlane.d
    end
  end
end