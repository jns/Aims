
module Aims
  
  # Utility class for parsing an Aims geometry file
  # Example Usage:
  #  uc = Aims::GeometryParser.parse("geometry.in")
  class GeometryParser
    
    # Parse a String representation of a geometry.in file
    # - +str+ The String to parse
    # - Return the Aims::Geometry object that was parsed
    def GeometryParser.parse_string(str)
      GeometryParser.parse_io(str)
    end
    
    # Parse an IO object representation of a geometry.in file
    # - +io+ The IO object to parse
    # - Return the Aims::Geometry object that was parsed
    def GeometryParser.parse_io(io)
      atoms = Array.new
      vectors = nil
      io.each_line{|line|
        case line
        when /\w*#.*/
          # Comment line, Do nothing
        when /atom/
          a, x, y, z, species = line.split(' ')
		  atom = Atom.new(x.to_f,y.to_f,z.to_f,species)
          atoms << atom
        when /lattice_vector/
          a, x, y, z = line.split(' ')
          vectors = Array.new if vectors.nil?
          vectors << Vector[x.to_f,y.to_f,z.to_f]
        when /constrain_relaxation/
          a, c = line.split(' ')
          atoms.last.constrain << c
        end
      }
      Geometry.new(atoms, vectors)
    end

    # Parse a geometry.in file
    # - +filename+ the file to parse
    # - return the Aims::Geometry object
    def GeometryParser.parse(filename)
      f = File.open(filename, 'r')
      cell = GeometryParser.parse_io(f)
      f.close
      return cell
    end
    
  end
  
end


