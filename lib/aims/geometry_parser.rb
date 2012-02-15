
module Aims
  
  class GeometryParser
        
    def GeometryParser.parse_string(str)
      GeometryParser.parse_io(str)
    end
    
    def GeometryParser.parse_io(io)
      atoms = Array.new
      vectors = nil
      io.each_line{|line|
        case line
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
      UnitCell.new(atoms, vectors)
    end

    def GeometryParser.parse(filename)
      f = File.open(filename, 'r')
      cell = GeometryParser.parse_io(f)
      f.close
      return cell
    end
    
  end
  
end


