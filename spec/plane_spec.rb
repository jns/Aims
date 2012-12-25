require 'aims'
include Aims

describe Plane do
  p = Plane.new(0,0,1,0)
  
  it "should intersect the ray +z at (0 0 0)" do
    a = Vector[0,0,0]
    b = Vector[0,0,1]
    x = p.intersection_with_ray(a, b)
    x.should eq(a)
  end
  
  it "should intersect the ray +z with origin (0,0,-1) at (0 0 0)" do
    a = Vector[0,0,-1]
    b = Vector[0,0,1]
    x = p.intersection_with_ray(a, b)
    x.should eq(Vector[0,0,0])
  end
  
  it "should not intersect the ray +x" do
    a = Vector[0,0,0]
    b = Vector[1,0,0]
    x = p.intersection_with_ray(a,b)
    x.should be_nil
  end
  
  it "should not intersect the +z ray with origin (0,0,1)" do
    a = Vector[0,0,1]
    b = Vector[0,0,1]
    x = p.intersection_with_ray(a,b)  
    x.should be_nil  
  end
  
end