#!/usr/bin/env ruby

require 'aims'
files = ARGV
if files.empty?
  puts "usage: #{File.basename $0} file1 [file2] ..."
  exit
end

STDOUT.sync = true
puts "%-10s \t %-20s \t %-15s \t %9s \t %7s \t %10s \t %12s \t %8s \t %10s" % %w(RUN FILE TOTAL_ENERGY NUM_ATOMS K-GRID CONVERGED RELAX_STEPS SC_ITERS TOTAL_TIME)
format = "%-10s \t %-20s \t %+15f \t %9i \t %7s \t %10s \t %12i \t %8i \t %10.2f"
files.each{|f|
  run = f.split(".").last
  begin
    o = Aims::OutputParser.parse(f)
    puts format % [run, f[0...20], 
                   (o.total_energy.nan? ? Float::NAN : o.total_energy), 
                   (o.n_atoms or -1), 
                   (o.k_grid ? o.k_grid.squeeze : "-"), 
                   o.geometry_converged, 
                   (o.n_relaxation_steps or -1), 
                   (o.n_sc_iterations or -1),
                   o.total_cpu_time]
  rescue 
    puts [run, f, "***ERROR***", $!.message].join("\t")
  end
}
