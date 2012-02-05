#!/usr/bin/ruby

require "Aims.rb"
require 'mathn'

module Aims
  
  class GeometryStep
    attr_reader :sc_iterations
    attr_accessor :geometry, :total_energy, :total_corrected_energy, :chemical_potential, :forces

    def initialize
      @sc_iterations = Array.new
      @forces = Array.new
    end

    def sc_iteration
      self.sc_iterations.last
    end

    def total_energy_per_atom
      self.total_energy/geometry.size rescue "N/A"
    end
    
    def total_corrected_energy_per_atom
      self.total_corrected_energy/geometry.size rescue "N/A"
    end
    
    def total_cpu_time
      val = self.sc_iterations.inject(0){|total, iter|
        total = total + iter.total_cpu_time
      }
      val or 0
    end
    
    def total_wall_time
      val = self.sc_iterations.inject(0){|total, iter|
        total = total + iter.total_wall_time
      }
      val or 0 
    end
  end
  
  class SCIteration
    attr_accessor :d_etot, :d_eev, :d_rho, :timings

    def initialize
      self.timings = Array.new
    end
    
    def total_cpu_time
      begin
        total= self.timings.find{|t| t[:description] =~ /Time for this iteration/}
        total[:cpu_time]
      rescue
        self.timings.inject(0){|total, time|
          total = total + (time[:cpu_time] || 0)
        }
      end
    end

    def total_wall_time
      begin
        time = self.timings.find{|t| t[:description] =~ /Time for this iteration/}
        time[:wall_time] or 0
      rescue
        self.timings.inject(0){|total, time|
          total = total + (time[:wall_clock_time] || 0)
        }
      end
    end
  end
  
  class AimsOutput
    attr_accessor :geometry_steps, :k_grid, :original_file, :geometry_converged, :timings, :computational_steps, :n_atoms
    
    def initialize
      self.geometry_steps = Array.new
      self.geometry_converged = false
      self.timings = {}
    end
    
    # Returns the best available value of the total energy
    def total_energy
      etot = self.geometry_steps.collect{|gs| gs.total_energy }.compact.last
      if etot.nil? 
        "NaN"
      else
        etot
      end
    end
    
    def final_geometry
      self.geometry_steps.last.geometry
    end
    
    def final_step
      self.geometry_steps.last
    end
    
    def geometry_step
      self.geometry_steps.last
    end
    
    def n_relaxation_steps
      self.geometry_steps.size - 1
    end
    
    def n_sc_iterations
      self.geometry_steps.inject(0){|total, step|
        total = total + step.sc_iterations.size
      }
    end

    def total_cpu_time
      begin 
        self.timings.find{|t| t[:description] =~ /Total time$/}[:cpu_time]
      rescue
        self.geometry_steps.inject(0) {|total, step|
          total = total + step.total_cpu_time
        }
      end
    end
    
    def sc_iteration
      self.geometry_step.sc_iteration
    end
  end
  
  class OutputParser
    
    def OutputParser.parse_input_geometry(io, n_atoms)
      atoms = []
      n_atoms.times do 
        fields = io.readline.split(' ')
        a = Atom.new
        a.x, a.y, a.z = fields[4].to_f, fields[5].to_f, fields[6].to_f
        a.species = fields[3]
        atoms << a
      end
      UnitCell.new(atoms)
    end

    def OutputParser.parse_updated_geometry(io, n_atoms)
      vectors = []
      3.times do 
        fields = io.readline.split(' ')
        vectors << [fields[1].to_f, fields[2].to_f, fields[3].to_f]
      end
      io.readline
      atoms = []
      n_atoms.times do 
        fields = io.readline.split(' ')
        a = Atom.new
        a.x, a.y, a.z = fields[1].to_f, fields[2].to_f, fields[3].to_f
        a.species = fields[4]
        atoms << a
      end
      UnitCell.new(atoms, vectors)
    end
    
    def OutputParser.parse_sc_timings(io)
      line = io.readline
      timings = []
      until line =~ /---/
        desc, times = line.split(":")
        fields = times.split(" ")
        cpu = fields[0].to_f
        wall = fields[2].to_f
        timings << {:description => desc.strip, :cpu_time => cpu, :wall_clock_time => wall}
        line = io.readline
      end
      timings
    end
    
    def OutputParser.parse_detailed_time_accounting(io)
      line = io.readline
      timings = []
      desc, times = line.split(":")
      until times.nil?
        fields = times.split(" ")
        cpu = fields[0].to_f
        wall = fields[2].to_f
        timings << {:description => desc.strip[1..-1].strip, :cpu_time => cpu, :wall_clock_time => wall}
        line = io.readline
        desc, times = line.split(":")
      end
      timings 
    end
    
    def OutputParser.parse_computational_steps(io)
      line = io.readline
      steps = []
      desc, value = line.split(":")
      until value.nil?
        steps << {:description => desc.strip[1..-1].strip, :value => value.to_f}
        line = io.readline
        desc, value = line.split(":")
      end
      steps 
    end

    def OutputParser.parse(filename)
  
      n_atoms = 0
      vectors = []
      retval = AimsOutput.new
      retval.original_file = filename
      
      File.open(filename, 'r') do |f|
        f.each_line{|line|
          case line
            when /Found k-point grid:/
              retval.k_grid = line.split(":")[1].strip
              
            when /Computational steps:/
              retval.computational_steps = OutputParser.parse_computational_steps(f)
            
            when /Detailed time accounting/
              retval.timings = OutputParser.parse_detailed_time_accounting(f)
              
            when /Begin self-consistency iteration/
              retval.geometry_step.sc_iterations << SCIteration.new

            when /End self-consistency iteration/
              retval.sc_iteration.timings = OutputParser.parse_sc_timings(f)

            when /Change of charge density/
              retval.sc_iteration.d_rho = line.split(' ')[6].to_f
          
            when /Change of sum of eigenvalues/
              retval.sc_iteration.d_eev = line.split(' ')[7].to_f
              
            when /Change of total energy/
              retval.sc_iteration.d_etot = line.split(' ')[6].to_f
            
            when /\|\ Total energy corrected/
              retval.geometry_step.total_corrected_energy = line.split(' ')[5].to_f
            
            when /\|\ Total energy uncorrected/
              retval.geometry_step.total_energy = line.split(' ')[5].to_f
              
            when /\|\ Number\ of\ atoms/
                n_atoms = line.split(' ')[5].to_i
                retval.n_atoms = n_atoms
                
            when  /\|\ Chemical potential/
              retval.geometry_step.chemical_potential = line.split(' ')[8].to_f
                
            when /Input\ geometry\:/
              line = f.readline
              if line=~/\|\ Unit\ cell\:/
                  3.times { 
                    line = f.readline
                    fields = line.split(' ')
                    vectors << Vector[fields[1].to_f, fields[2].to_f, fields[3].to_f]
                  }
              end
              2.times {f.readline} 
              retval.geometry_steps << GeometryStep.new 
              retval.geometry_step.geometry = OutputParser.parse_input_geometry(f, n_atoms)
              retval.geometry_step.geometry.lattice_vectors = vectors

            when /\ Updated\ atomic\ structure\:/
              f.readline
              retval.geometry_steps << GeometryStep.new 
              retval.geometry_step.geometry = OutputParser.parse_updated_geometry(f, n_atoms)
            #  retval.geometry_step.geometry.lattice_vectors = vectors
              
            when /\ Final\ atomic\ structure\:/
              f.readline
              retval.geometry_step.geometry = OutputParser.parse_updated_geometry(f, n_atoms)
             # retval.geometry_step.geometry.lattice_vectors = vectors
            when /\  Total\ atomic\ forces/
              line = f.readline
              until line =~ /---/
                  fields = line.split(' ')
                  retval.geometry_step.forces << Vector[fields[2].to_f, fields[3].to_f, fields[4].to_f]
                  line = f.readline
              end
            when /Present geometry is converged./
              retval.geometry_converged = true
          end
        }
      end

      return retval
      
    end
  end
  
end

