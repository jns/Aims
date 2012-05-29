require 'aims'
include Aims

describe OutputParser do
  
  updated_geometry_non_periodic =<<-EOL
                         x [A]             y [A]             z [A]
            atom         0.00000000       -0.03607483        0.00000000  O
            atom         0.74134141       -0.68896259        0.00000000  H
            atom        -0.74134141       -0.68896259        0.00000000  H
  
  
EOL
  
  updated_geometry_periodic =<<-EOL
                         x [A]             y [A]             z [A]
  lattice_vector         1.99266923        1.99264227       -0.00253414
  lattice_vector        -0.00256109        1.99264197        1.99264198
  lattice_vector         1.99266921       -0.00253413        1.99264226

            atom         0.00000000        0.00000000        0.00000000  Al

  Fractional coordinates:
                         L1                L2                L3
       atom_frac         0.00000000        0.00000000        0.00000000  Al  
EOL
 
 it "should have 3 atoms" do
   g = OutputParser.parse_updated_geometry(StringIO.new(updated_geometry_non_periodic), 3)
   g.atoms[0].should eq(Atom.new(0,-0.3607483, 0, "O"))
   g.atoms[1].should eq(Atom.new(0.74134141, -0.68896259, 0, "H"))
   g.atoms[2].should eq(Atom.new(-0.74134141, -0.68896259, 0, "H"))
 end
 
 it "should have 3 lattice vectors and 1 atom" do
   g = OutputParser.parse_updated_geometry(StringIO.new(updated_geometry_periodic), 3)
   g.lattice_vectors[0].should eq(Vector[1.99266923, 1.99264227, -0.00253414])
   g.lattice_vectors[1].should eq(Vector[-0.00256109, 1.99264197,  1.99264198])
   g.lattice_vectors[2].should eq(Vector[ 1.99266921, -0.00253413, 1.99264226])
   g.atoms[0].should eq(Atom.new(0,0,0,"Al"))
 end
 
end