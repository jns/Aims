require 'aims'
require "test/unit"

class TestVolume < Test::Unit::TestCase
  
  include Aims
  
  def setup
    # A unit cube
    @p1 = Plane.new(1,0,0,1,1,1)
    @p2 = Plane.new(0,1,0,1,1,1)
    @p3 = Plane.new(0,0,1,1,1,1)
    @p4 = Plane.new(-1,0,0,0,0,0)
    @p5 = Plane.new(0,-1,0,0,0,0)
    @p6 = Plane.new(0,0,-1,0,0,0)
    
    # A diagonal to make a wedge
    @p7 = Plane.new(1,1,0,1,0,0)
    
    # A plane just outside the unit cube
    @p8 = Plane.new(1,0,0,1.1, 0, 0)
  end
  
  def teardown
    
  end
  
  def test_intersections
    assert_equal(1, Volume.intersection_points([@p1, @p2, @p3]).size)
    assert_equal(8, Volume.intersection_points([@p1, @p2, @p3, @p4, @p5, @p6]).size)
    assert_equal(6, Volume.intersection_points([@p3, @p4, @p5, @p6, @p7]).size)
    assert_equal(8, Volume.intersection_points([@p1, @p2, @p3, @p4, @p5, @p6, @p8]).size)
  end
  
  def test_wedge
    wedge = Volume.new([@p3,@p4,@p5,@p6,@p7])
    assert(wedge.contains_point(0,0,0))
    assert((not wedge.contains_point(1,1,1)))
  end
  
  def test_cube
    cube = Volume.new([@p1, @p2, @p3, @p4, @p5, @p6, @p8])
    assert(cube.contains_point(0.5,0.5,0.5))
  end
  
  def test_bounding_box
    cube = Volume.new([@p1, @p2, @p3, @p4, @p5, @p6, @p8])
    bbox = cube.bounding_box
    assert(bbox.contains_point(1,1,1))
    assert(bbox.contains_point(0,0,0))
  end
  
  def test_fill
    s = 10
    cube = Volume.new([Plane.new(1,0,0,s,s,s),
                       Plane.new(0,1,0,s,s,s),
                       Plane.new(0,0,1,s,s,s),
                       Plane.new(-1,0,0,0,0,0),
                       Plane.new(0,-1,0,0,0,0),
                       Plane.new(0,0,-1,0,0,0)])
    
    zb = ZincBlende.new("Ga", "As", 5.75)
    assert_equal(1, zb.fill_volume(cube).atoms.size)
    
  end
end