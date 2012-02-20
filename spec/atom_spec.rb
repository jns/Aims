require 'aims'
include Aims

describe Atom do 

  a = Atom.new(0,0,0,"Atom")
  
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
  
end