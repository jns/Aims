module Aims
  
# A volume is defined by a minimum of four planes that intersect in a minimum 
# of four points. The normals of the planes must all point outward. This is tested
# 
  class Volume
    
    # Add Vectorize as class methods
    class<<self
      include Vectorize
    end
    
    # Quick recursive method for calculating combinations
    def Volume.choose(list, num, head = [])

      _head = head.dup
      _list = list.dup
      _num = num

      if _num == 0 
        return [_head]
      end
        
      new_heads = []
      while _list.size > _num-1
          h = _head + [_list.shift]
          new_heads += Volume.choose(_list, num-1, h)
        end
      return new_heads
      
    end
    
    
    # Return an array of tuples that define the 
    # vertices of intersection of these planes
    # Vertices are removed that lie in front of any plane
    # @param [Array<Plane>] planes An array of Planes that define the boundaries of the volume
    # @return [Array] An array of tuples that define the vertices of intersection of these planes
    def Volume.intersection_points(planes)
      combos = Volume.choose(planes, 3)
      points = []
      combos.each{|c|
        n1 = c[0].unit_normal
        n2 = c[1].unit_normal
        n3 = c[2].unit_normal
        d = Matrix[n1, n2,n3].transpose.det
        
        # The determinant is zero if any two planes are parallel
        unless (d == 0)
          p1 = c[0].any_point_on_plane
          p2 = c[1].any_point_on_plane
          p3 = c[2].any_point_on_plane
        
          # This defines the point of intersection of three planes.
          points << (cross(n2,n3)*dot(p1, n1) + cross(n3,n1)*dot(p2,n2) + cross(n1,n2)*dot(p3, n3))*(1/d)
        end
      }
      
      # Only keep the points that are behind all planes
      keepers = []
      points.each{|pt|
        keep = true
        planes.each {|pl|
          keep = (pl.distance_to_point(pt[0], pt[1], pt[2]) <= 0)
          break unless keep
        }
        keepers << pt if keep
      }
      return keepers
    end
    
    def initialize(planes)
      
      points = Volume.intersection_points(planes)
      if (4 > points.size)
        raise "Planes do not intersect in a closed volume."
      end
      @points = points
      @planes = planes
    end
    
    # Return the bounding box for this volume
    def bounding_box
      unless @bbox
        p = @points[0]
        minX = p[0]
        maxX = p[0]
        minY = p[1]
        maxY = p[1]
        minZ = p[2]
        maxZ = p[2]
        @points.each{|p|
          minX = p[0] if p[0] < minX  
          maxX = p[0] if p[0] > maxX  
          minY = p[1] if p[1] < minY  
          maxY = p[1] if p[1] > maxY  
          minZ = p[2] if p[2] < minZ  
          maxZ = p[2] if p[2] > maxZ  
        }
        @max = Vector[maxX, maxY,maxZ]
        @min = Vector[minX, minY, minZ]
        @bbox = Volume.new([Plane.new(-1,0,0, minX, minY, minZ), 
                           Plane.new(0,-1,0, minX, minY, minZ),
                           Plane.new(0,0,-1, minX, minY, minZ),
                           Plane.new(1,0,0, maxX, maxY, maxZ),
                           Plane.new(0,1,0, maxX, maxY, maxZ), 
                           Plane.new(0,0,1, maxX, maxY, maxZ)])
      end
      @bbox
    end
    
    # Return the point on the bounding box that represents the maximum value
    # of any cartesian coordinate (the upper right corner)
    def max_point
      # generate the bounding box if not already done
      bounding_box
      # return the max
      @max
    end
    
    # Return the point on the bounding box that represents the minimum value
    # of any cartesian coordinate (the lower left corner).
    def min_point
      # generate the bounding box if not already done
      bounding_box
      # return the min
      @min
    end
    
    # A volume contains a point if it lies behind all the planes
    def contains_point(x,y,z)
      behind = true
      @planes.each{|p|
        behind = (0 >= p.distance_to_point(x,y,z))
        break if not behind
      }
      return behind
    end
    
  end

end