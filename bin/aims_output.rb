#!/usr/bin/env ruby

require 'optparse'
require 'aims'

options = {:step => :all}

OptionParser.new do |opts|
  opts.on('-s', '--step [N]', 'Output information for relaxation step.', 
                                  "Specify an integer, 'first', 'last', or 'all'", 
                                  "Default is 'all'") do |s|
    case s
    when /([1-9]+)/
      options[:step] = $1.to_i
    when "first"
      options[:step] = :first
    when "last"
      options[:step] = :last
    else
      options[:step] = :all
    end
  end

  opts.on('-d', '--geometry-delta', 'Output change from geometry in to geometry final') do
    options[:geometry_delta] = true
  end

  opts.on('-c', '--self-consistency', 'Output self-consistency information') do
    options[:self_consistency] = true
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
    
    steps = case options[:step]
    when Integer
      stepno = options[:step]
      if stepno < 0
        [output.geometry_steps.last]
      elsif stepno < output.geometry_steps.size
        [output.geometry_steps[stepno]]
      else
        [output.geometry_steps.last]
      end
    when :first
      [output.geometry_steps.first]
    when :last
      [output.geometry_steps.last]
    else
      output.geometry_steps
    end
    
    steps.each_with_index{|step, i|

      total_relaxations += 1
      total_sc_iterations += step.sc_iterations.size

      sciter_format  = "%-20s %20i"
      timings_format = "%-35s %20.5f"
      energy_format  = "%-35s %20.5e"
      
      puts "= Relaxation Step #{step.step_num} ="
      
      indent = "  "
      puts indent + sciter_format % ["SC Iterations", step.sc_iterations.size]
      puts indent + energy_format % ["Total Energy", step.total_energy]
      puts indent + timings_format % ["Total CPU time", step.total_cpu_time]        
      puts indent + timings_format % ["Total Wall time", step.total_wall_time]        
      unless step.forces.empty?
        puts "Max force: #{step.forces.max{|a,b| a.r <=> b.r}.r}" if options[:forces]
      end
      if options[:timings]
        puts "  Cumulative SC Timings:"
        step.timings.each{|t| puts "   " +timings_format % [t[:description], t[:cpu_time]]} 
      end
      
      if options[:self_consistency]
        
        indent = "    "
        
        # Iterate over each sc iteration
        step.sc_iterations.each_with_index{|sc_iter, iter|
          # SC Iteration Header
          puts "  == SC Iteration #{iter} =="

          # Output convergence criterion
          puts indent + energy_format % ["Change in total energy", sc_iter.d_etot]
          puts indent + energy_format % ["Change in sum of eigenvalues", sc_iter.d_eev]
          puts indent + energy_format % ["Change in charge density", sc_iter.d_rho]
          
          # Output timings if requested
          if options[:timings]
            if sc_iter.timings
              sc_iter.timings.each{|t|
                puts indent + timings_format % [t[:description], t[:cpu_time]]
              } 
            else
              puts "No timing data available."
            end
          end
          puts ""
        }
      end
      
      puts step.geometry.format_geometry_in if options[:geometry]

      puts "\n\n"
      
      
    }

    if options[:geometry_delta]
      puts "= Change in atomic positions for calculation"
      puts output.geometry_steps.last.geometry.delta(output.geometry_steps.first.geometry)
      unless output.geometry_converged
        puts "Warning Geometry not converged!"
      end
    end
  }



  # puts "Total relaxation steps: #{total_relaxations}"
  # puts "Total sc iterations: #{total_sc_iterations}"

rescue
  puts ""
  puts "Sorry. There was an error parsing the remainder of the file."
  puts ""
  exit
end

