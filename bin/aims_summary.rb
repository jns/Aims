#!/usr/bin/env ruby

require 'aims'
files = ARGV
STDOUT.sync = true
puts "%-10s | %-20s | %-15s | %9s | %7s | %10s | %12s | %8s | %10s" % %w(RUN FILE TOTAL_ENERGY NUM_ATOMS K-GRID CONVERGED RELAX_STEPS SC_ITERS TOTAL_TIME)
format = "%-10s | %-20s | %+15e | %9i | %7s | %10s | %12i | %8i | %10.2f"
files.each{|f|
  run = f.split(".").last
  begin
    o = Aims::OutputParser.parse(f)
    puts format % [run, f[0...20], (o.total_energy.nan? ? Float::NAN : o.total_energy), (o.n_atoms or -1), (o.k_grid ? o.k_grid.squeeze : "-"), o.geometry_converged, (o.n_relaxation_steps or -1), (o.n_sc_iterations or -1), o.total_cpu_time]
  rescue 
    puts [run, f, "***ERROR***", $!.message].join("\t")
  end
}
