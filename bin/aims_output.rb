#!/usr/bin/env ruby

require 'optparse'
require 'aims'



options = {:step => :all}

optParser = OptionParser.new do |opts|
  opts.banner = "usage: #{File.basename $0} [options] file1 [file2 ...]" 
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

  opts.on('--debug', 'Debug output') do 
    options[:debug] = true
  end

  opts.on('--geometry-delta', 'Display change from input geometry to final geometry') do
    options[:geometry_delta] = true
  end

  opts.on('-c', '--self-consistency', 'Output self-consistency information') do
    options[:self_consistency] = true
  end

  opts.on('-f', 'Output max force component for each geometry relaxation step') do 
    options[:forces] = true
  end

  opts.on('-t', 'Output timings') do
    options[:timings] = true
  end

end

begin
  optParser.parse!(ARGV)
  
  files = ARGV
  if files.empty?
    puts optParser.help
    exit
  end
  outputs = files.collect{|f|
      Aims::OutputParser.parse(f)
  }

  int_format = "%-20s %20i"
  float_format = "%-20s 20.5f"
  exp_format = "%-20s %20.5e"
  
  sciter_format  = "%-20s %20i"
  timings_format = "%-35s %20.5f"
  energy_format  = "%-35s %20.5f"
  force_format  = "%-35s %20.5e"

  total_sc_iterations = 0
  total_relaxations = 0

  outputs.each{|output|
    
    puts "**************************************************************************"
    puts "**"
    puts "**   #{output.original_file}"
    puts "**"
    puts "**************************************************************************"
    
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
      
      puts "= Relaxation Step #{step.step_num} ="
      
      indent = "  "
      puts indent + sciter_format % ["SC Iterations", step.sc_iterations.size]
      puts indent + energy_format % ["Total Energy", step.total_energy]
      puts indent + timings_format % ["Total CPU time", step.total_cpu_time]        
      puts indent + timings_format % ["Total Wall time", step.total_wall_time]        
      if options[:forces] and not step.forces.empty?
        puts indent + force_format % ["Max Force", step.forces.max{|a,b| a.r <=> b.r}.r] 
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
          puts indent + exp_format % ["Change in total energy", sc_iter.d_etot]
          puts indent + exp_format % ["Change in sum of eigenvalues", sc_iter.d_eev]
          puts indent + exp_format % ["Change in charge density", sc_iter.d_rho]
          
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

      puts "\n\n"
      
      
    }

    puts "= Calculation Summary ="
    unless output.geometry_converged
      puts "Warning Geometry not converged!"
    end

    puts int_format % ["Number of SC Iterations found:", total_sc_iterations]
    puts int_format % ["Number of Relaxation Steps found:", total_relaxations]
    puts float_format % ["Total CPU time", output.total_cpu_time]   

    output.computational_steps.each{|cs| 
      puts int_format % [cs[:description], cs[:value]]
    }
    if options[:timings]
      output.timings.each{|t| puts "   " +timings_format % [t[:description], t[:cpu_time]]} 
    end
    # puts timings_format % ["Total Wall time", output.total_wall_time]        

    if options[:geometry_delta]
      puts "= Change in atomic positions for calculation"
      puts output.geometry_steps.last.geometry.delta(output.geometry_steps.first.geometry)
    end
  }



  # puts "Total relaxation steps: #{total_relaxations}"
  # puts "Total sc iterations: #{total_sc_iterations}"

rescue
  puts ""
  puts "Sorry. There was an error parsing the remainder of the file."
  if options[:debug]
    puts $!.message 
    puts $!.backtrace
  else
    puts "Rerun with --debug for more info"
  end
  puts ""
  exit
end

