
# Copyright 2012 Joshua Shapiro
# released under the MIT license

require 'aims/atom.rb'
require 'aims/bond.rb'
require 'aims/geometry_parser.rb'
require 'aims/output.rb'
require 'aims/plane.rb'
require 'aims/unit_cell.rb'
require 'aims/volume.rb'
require 'aims/wurtzite.rb'
require 'aims/zinc_blende.rb'

module Aims
  # Define utility methods at the module level
  
  # Dot product of two n-element arrays
  def dot(a, b)
    unless a.size == b.size
      raise "Vectors must be the same length"
    end
    
    # Make element-by-element array of pairs
    (a.to_a).zip(b.to_a).inject(0) {|tot, pair| tot = tot + pair[0]*pair[1]}
  end

  # Cross product of two arrays of length 3
  def cross(b,c)
    unless b.size == 3 and c.size == 3
      raise "Vectors must be of length 3"
    end
    Vector[b[1]*c[2] - b[2]*c[1], b[2]*c[0] - b[0]*c[2], b[0]*c[1] - b[1]*c[0]]
  end
end