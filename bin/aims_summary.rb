#!/usr/bin/env ruby

require 'aims'
files = ARGV
STDOUT.sync = true
puts ["RUN", "FILE", "TOTAL_ENERGY", "NUM_ATOMS", "K-GRID", "CONVERGED", "RELAXATION_STEPS", "SC_ITERATIONS","TOTAL_TIME"].join("\t")
files.each{|f|
  run = f.split(".").last
  begin
    o = Aims::OutputParser.parse(f)
    puts [run, f, o.total_energy, o.n_atoms, (o.k_grid or "-"), o.geometry_converged, o.n_relaxation_steps, o.n_sc_iterations, o.total_cpu_time].join("\t")
  rescue 
    puts [run, f, "***ERROR***", $!.to_s].join("\t")
	puts $!.backtrace
  end
}
