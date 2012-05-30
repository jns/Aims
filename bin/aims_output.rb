#!/usr/bin/env ruby

require 'optparse'
require 'aims'

options = {}
OptionParser.new do |opts|
  opts.on('-g', '--geometry', 'Output geometry information') do
    options[:geometry] = true
  end

  opts.on('-d', '--geometry-delta', 'Output change from geometry in to geometry final') do
    options[:geometry_delta] = true
  end

  opts.on('-G [list of iterations]', '--geometry-compare', Array, 'Compare converged geometry for two files') do |list|
    options[:geometry_compare] = list
  end
  
  opts.on('-s', '--self-consistency', 'Output self-consistency information') do
    options[:self_consistency] = true
  end
  
  opts.on('-e', 'Output energy information') do 
    options[:energy] = true
  end

  opts.on('-E n', Integer, 'Output total energy for given iteration') do |n|
     options[:energy] = n
  end  
  opts.on('-f', 'Output forces') do 
    options[:forces] = true
  end
  
  opts.on('-t', 'Output timings') do
    options[:timings] = true
  end
  
end.parse!

begin
  
  files = ARGV
  outputs = files.collect{|f|
      Aims::OutputParser.parse(f)
  }
 
  total_sc_iterations = 0
  total_relaxations = 0
  
    outputs.each{|output|
      puts output.original_file
      output.geometry_steps.each_with_index{|step, i|
        
        total_relaxations += 1
        total_sc_iterations += step.sc_iterations.size
       
        puts "Total Energy: #{step.total_energy}" if options[:energy] == i 
        puts "Iteration #{i}" if options[:geometry] or options[:energy] == true
        puts "Total Corrected Energy #{step.total_corrected_energy}" if options[:energy] == true
        puts "Total Corrected Energy/atom #{step.total_corrected_energy_per_atom}" if options[:energy] == true
        puts "Chemical Potential #{step.chemical_potential}" if options[:energy] == true
        puts "CPU Time: #{step.total_cpu_time},  Wall Time: #{step.total_wall_time}" if options[:timings]
        unless step.forces.empty?
          puts "Max force: #{step.forces.max{|a,b| a.r <=> b.r}.r}" if options[:forces]
        end
        puts step.geometry.format_geometry_in if options[:geometry]
      }
      
      if options[:geometry_delta]
        puts "Change"
        puts output.geometry_steps.last.geometry.delta(output.geometry_steps.first.geometry)
        unless output.geometry_converged
          puts "Warning Geometry not converged!"
        end
      end
    }
 
 if options[:geometry_compare] && outputs.size == 2
   list = options[:geometry_compare]
   puts list[0].to_i || "NA"
   puts list[1].to_i || "NA"
   g1 = outputs[0].geometry_steps[list[0].to_i].geometry || outputs[0].final_geometry
   g2 = outputs[1].geometry_steps[list[1].to_i].geometry || outputs[1].final_geometry
   puts "----#{outputs[0].original_file} vs. #{outputs[1].original_file}----"
   puts g2.delta(g1).format_geometry_in
 end
 
 if options[:self_consistency]
   outputs.each{|output|
     puts output.original_file
    output.geometry_steps.each_with_index{|step,i|
      puts "RelaxationStep\t#{i}"
      puts "Total Corrected Energy\t#{step.total_corrected_energy}"
      puts "Total Corrected Energy/atom\t#{step.total_corrected_energy_per_atom}"

      format = "%4i % 6e % 6e % 6e"
      puts "%4s %13s %13s %13s" % %w(iter d_etot d_eev d_rho)
      step.sc_iterations.each_with_index{|sc_iter, iter|
        puts format % [iter, sc_iter.d_etot, sc_iter.d_eev, sc_iter.d_rho]
        if options[:timings]
          sc_iter.timings.each{|t|
            puts "#{t[:description]}\t#{t[:cpu_time]}"
          }
        end
      }
    } 
  }
  
  if options[:timings]
    
  end
 end

 
 # puts "Total relaxation steps: #{total_relaxations}"
 # puts "Total sc iterations: #{total_sc_iterations}"

 rescue
   puts $!.message
   exit
 end
 
