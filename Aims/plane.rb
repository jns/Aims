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
      @d = a*x + b*y + c*z
    end

    # The distance to a point (x,y,z) is
    # D - Ax - By - Cz
    def distance_to_point(x, y, z)
      a*x + b*y + c*z - d
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