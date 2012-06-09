
require 'aims'
include Aims

describe ZincBlende do 

  zb = ZincBlende.new("Cation", "Anion", 5.75)

  context "Bulk Geometry" do 
    bulk = zb.get_bulk
    it "should have 2 atoms" do
      bulk.atoms.size.should eq(2)
    end
  end
  
  context "001 Slab" do 
    (1..6).each do |i|
      context "#{i} layers" do
        s = zb.get_001_surface(i, 20)
        it "should have #{i} atoms" do 
          s.atoms.size.should eq(i)
        end
    
        it "should have 0 constrained atoms" do 
          s.atoms.find_all{|a| a.constrained?}.size.should eq(0)
        end
      end
    end

    6.times do |i|
      context "#{i} constrained layers" do 
        s = zb.get_001_surface(6, 20, i)
        it "should have 6 atoms" do 
          s.atoms.size.should eq(6)
        end
    
        it "should have 20AA of vacuum" do
          maxz = s.atoms.max{|a,b| a.z <=> b.z}.z
          minz = s.atoms.min{|a,b| a.z <=> b.z}.z
          z_vector = s.lattice_vectors[2][2]
          (z_vector - (maxz - minz)).should eq (20)
        end
    
        it "should have i constrained atoms" do 
          s.atoms.find_all{|a| a.constrained?}.size.should eq(i)
        end
      end
    end
  end

  context "110 Slab" do 
    
    (1..5).each do |i|
      s = zb.get_110_surface(i, 20)
      it "should have #{i*2} atoms" do
        # The 110 surface has 2 atoms per monolayer
        s.atoms.size.should eq(i*2)
      end

      it "should have 20AA of vacuum" do
        maxz = s.atoms.max{|a,b| a.z <=> b.z}.z
        minz = s.atoms.min{|a,b| a.z <=> b.z}.z
        z_vector = s.lattice_vectors[2][2]
        (z_vector - (maxz - minz)).should eq (20)
      end
    end

    5.times do |i|
      it "Should have #{i*2} constrained atoms" do
        s = zb.get_110_surface(5, 20, i)
        s.atoms.find_all{|a| a.constrained? }.size.should eq(i*2)
      end
    end
  end

  context "111 Slab" do 
    s = zb.get_110_surface(5, 20)
    it "should have 10 atoms" do
      # The 110 surface has 2 atoms per monolayer
      s.atoms.size.should eq(10)
    end
    it "should have 20AA of vacuum" do
      maxz = s.atoms.max{|a,b| a.z <=> b.z}.z
      minz = s.atoms.min{|a,b| a.z <=> b.z}.z
      z_vector = s.lattice_vectors[2][2]
      (z_vector - (maxz - minz)).should eq (20)
    end
  end

end