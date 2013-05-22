Gem::Specification.new do |s|
  s.author = "Joshua Shapiro"
  s.email = "joshua.shapiro@gmail.com"
  s.description = "Support for generation and parsing of input and output files for FHI-AIMS DFT package"
  s.files = Dir.glob("{bin,lib}/**/*.rb") + %w(README.rdoc)
  s.homepage = "https://github.com/jns/Aims"
  s.name = "aims"
  s.require_path = 'lib'
  s.summary =<<-EOF
  This gem offers support for parsing and generating geometry and control files, 
  and parsing output files for the FHI-AIMS DFT code.
EOF
  s.version = "0.3.0"
  s.executables = ["aims_output.rb", "aims_summary.rb"]
  s.test_files = Dir.glob("spec/*.rb")
  s.add_development_dependency 'rspec'
end
