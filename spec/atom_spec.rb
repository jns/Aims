require 'aims'
include Aims

describe Atom do 

  a = Atom.new(0,0,0,"Atom")
  a_copy = a.copy
  
  it "Should not be constrained" do 
    a.constrained?.should be_false
  end
  
  it "constrained should be TRUE" do 
    a.constrain = true
    a.constrained?.should be_true
  end

  it "constrained should be '.true.'" do 
    a.constrain = ".true."
    a.constrained?.should be_true
  end
  
  it "should == its copy" do
    (a == a_copy).should be_true
  end
  
  it "should have the same hash as its copy" do
    a.hash.should eq(a_copy.hash)
  end
  
  it "should eql? its copy" do
    (a.eql?(a_copy)).should be_true
  end
  
  context "An array of duplicates" do
    it "should have one element when uniq" do
      [a, a_copy].uniq.size.should eq(1)
    end
  end
end