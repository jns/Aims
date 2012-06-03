module Aims

  # A geometry relaxation step in the calculation
  class GeometryStep

    # The relaxation step number
    attr_accessor :step_num

    # An Array of Aims::SCIteration 
    attr_reader :sc_iterations

    # The Aims::Geometry
    attr_accessor :geometry

    # The total energy for this geometry as a Float
    attr_accessor :total_energy

    # The total corrected energy for this geometry as a Float
    attr_accessor :total_corrected_energy
    
    # The chemical potential for this geometry
    attr_accessor :chemical_potential
    
    # The forces for this geometry
    attr_accessor :forces

    # Initialize a new geometry step
    def initialize
      @sc_iterations = Array.new
      @forces = Array.new
      @step_num = -1
      @total_energy = Float::NAN
      @total_corrected_energy = Float::NAN
      @chemical_potential = Float::NAN
    end

    # Return the last self-consistency iteration for this relaxation step
    def sc_iteration
      self.sc_iterations.last
    end

    # Return the total energy per atom for this relaxation step
    def total_energy_per_atom
      self.total_energy/geometry.size rescue "N/A"
    end
    
    # Return the total corrected energy per atom for this relaxation step
    def total_corrected_energy_per_atom
      self.total_corrected_energy/geometry.size rescue "N/A"
    end
    
    # A hash with keys :description, :cpu_time, and :wall_time
    # with detailed time accounting for this geometry step.
    def timings
      timings = Timings.new
      self.sc_iterations.each{|sc_iter|
        timings.add!(sc_iter.timings)
      }
      timings
    end

    # The total CPU time deterimed by summing the time for
    # each self-consistency iteration
    def total_cpu_time
      self.sc_iterations.inject(0){|tot, sc_iter|
        tot += sc_iter.total_cpu_time
      }
    end
    
    
    # The total wall clock time deterimed by summing the time for
    # each self-consistency iteration
    def total_wall_time
      self.sc_iterations.inject(0){|tot, sc_iter|
        tot += sc_iter.total_wall_time
      }
    end
  end
  
  # A class for encapsulating computational timing information
  class Timings
    
    # Initialize a new timing 
    def initialize
      @timing_hash = {}
    end
    
    # Enumerate over an array of hashes that looks like:
    # [{:description => "Something", :cpu_time => 10.0, :wall_time => 10.1}, ...]
    def each
      @timing_hash.each_pair{|desc,timings|
        h = {:description => desc, :cpu_time => timings[:cpu_time], :wall_time => timings[:wall_time]}
        yield h
      }
    end
    
    # Add another timings object to this one
    def add!(timings)
      timings.descriptions.each{|d|
        add_cpu_time(d, timings.cpu_time(d))
        add_wall_time(d, timings.wall_time(d))
      }
    end
    
    # Get all the descriptions for this timing object
    def descriptions
      @timing_hash.keys
    end
    
    # Add  cpu timing data for the given description
    def add_cpu_time(desc, time)
      @timing_hash[desc] = {:cpu_time => 0, :wall_time => 0} unless @timing_hash[desc]
      @timing_hash[desc][:cpu_time] += time
    end
    
    # Add  wall timing data for the given description
    def add_wall_time(desc, time)
      @timing_hash[desc] = {:cpu_time => 0, :wall_time => 0} unless @timing_hash[desc]
      @timing_hash[desc][:wall_time] += time
    end
    
    # Get wall timing data for the given description
    def wall_time(desc)
      if @timing_hash[desc]
        @timing_hash[desc][:wall_time] || 0
      else
        0
      end
    end
    
    # Get the total wall time
    def total_wall_time
      @timing_hash.inject(0) {|total, a|
        total += a[1][:wall_time]
      }
    end
    
    # Get cpu timing data for the given description
    def cpu_time(desc)
      if @timing_hash[desc]
        @timing_hash[desc][:cpu_time] || 0
      else
        0
      end
    end

    # Get the total cpu time
    def total_cpu_time
      @timing_hash.inject(0) {|total, a|
        total += a[1][:cpu_time]
      }
    end
  end
  
  # A single self-consistency iteration
  class SCIteration
    
    # The change in total energy
    attr_accessor :d_etot
    
    # The change in the sum of eigenvalues
    attr_accessor :d_eev
    
    # The change in charge density
    attr_accessor :d_rho
    
    # The Timings data for this iteration
    attr_accessor :timings

    # Initialize a new self-consistency iteration
    def initialize
      self.timings = Timings.new
      self.d_eev = 0
      self.d_rho = 0
      self.d_etot = 0
    end
    
    # Return the total CPU time for this iteration
    def total_cpu_time
      self.timings.cpu_time("Time for this iteration")
      # begin
      #   total= self.timings.find{|t| t[:description] =~ /Time for this iteration/}
      #   total[:cpu_time]
      # rescue
      #   self.timings.inject(0){|total, time|
      #     total = total + (time[:cpu_time] || 0)
      #   }
      # end
    end

    # Return the total wall clock time for this iteration
    def total_wall_time
      self.timings.wall_time("Time for this iteration")
      # begin
      #   time = self.timings.find{|t| t[:description] =~ /Time for this iteration/}
      #   time[:wall_time] or 0
      # rescue
      #   self.timings.inject(0){|total, time|
      #     total = total + (time[:wall_clock_time] || 0)
      #   }
      # end
    end
  end
  
  # An object encapsulating the data that is output from AIMS.
  # This object is generated from Aims::OutputParser.parse(filename)
  class AimsOutput
    # Each Aims::GeometryStep geometry relaxation step
    attr_accessor :geometry_steps
    
    # The k-point grid for periodic calculations
    attr_accessor :k_grid
    
    # The name of the calculation output file
    attr_accessor :original_file
    
    # Boolean, true if the geometry is converged
    attr_accessor :geometry_converged
    
    # The detailed time accounting data as a Aims::Timings
    # will be nil if the calculation did not complete
    attr_accessor :timings
    
    # ?
    attr_accessor :computational_steps
    
    # The number of atoms in the computation
    attr_accessor :n_atoms
    
    def initialize
      self.geometry_steps = Array.new
      self.geometry_converged = false
      self.timings = {}
    end
    
    # Returns the best available value of the total energy
    def total_energy
      etot = self.geometry_steps.collect{|gs| gs.total_energy }.compact.last
      if etot.nil? 
        Float::NAN
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
  
  # Parse an AIMS output file and generate an AimsOutput object
  # Invoke with 
  #  output = Aims::OutputParser.parse(filename)
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
      Geometry.new(atoms, nil, :dont_make_bonds)
    end

    def OutputParser.parse_atom_frac(line)
      nil
    end
    
    def OutputParser.parse_atom(line)
      fields = line.split(' ')
      a = Atom.new
      a.x, a.y, a.z = fields[1].to_f, fields[2].to_f, fields[3].to_f
      a.species = fields[4]
      a
    end

    def OutputParser.parse_lattice_vector(line)
      fields = line.split(' ')
      [fields[1].to_f, fields[2].to_f, fields[3].to_f]      
    end

    def OutputParser.parse_updated_geometry(io, n_atoms)
      io.readline
      
      vectors = []
      atoms = []
      continue = TRUE
      
      begin
        line = io.readline
        case line
        when /lattice_vector/
          vectors << OutputParser.parse_lattice_vector(line)
        when /atom\b/
          atoms << OutputParser.parse_atom(line)
        when /atom_frac/
          OutputParser.parse_atom_frac(line)
        when /^\s*$/
          # do nothing
        else 
          continue = FALSE
        end
      end while continue 
      
      vectors = nil if vectors.empty?
      
      Geometry.new(atoms, vectors, :dont_make_bonds)
    end
    
    def OutputParser.parse_sc_timings(io)
      line = io.readline
      timings = Timings.new
      until line =~ /---/
        desc, times = line.split(":")
        fields = times.split(" ")
        description = desc.sub("|", "").strip
        cpu = fields[0].to_f
        wall = fields[2].to_f
        timings.add_cpu_time(description, cpu)
        timings.add_wall_time(description, wall)
        line = io.readline
      end
      timings
    end
    
    def OutputParser.parse_detailed_time_accounting(io)
      line = io.readline
      timings = Timings.new
      desc, times = line.split(":")
      until times.nil?
        fields = times.split(" ")
        description = desc.sub("|", "").strip
        cpu = fields[0].to_f
        wall = fields[2].to_f
        timings.add_cpu_time(description, cpu)
        timings.add_wall_time(description, wall)
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
        steps << {:description => desc.sub("|", " ").strip, :value => value.to_f}
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
              
            when /Begin self-consistency loop/
              retval.geometry_step.sc_iterations << SCIteration.new
              
            when /Begin self-consistency iteration/
              retval.geometry_step.sc_iterations << SCIteration.new
                            
            when /End self-consistency iteration/, /End scf initialization - timings/
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
              retval.geometry_step.step_num = 0
              retval.geometry_step.geometry = OutputParser.parse_input_geometry(f, n_atoms)
              retval.geometry_step.geometry.lattice_vectors = vectors

            when /\ Updated\ atomic\ structure\:/
              last_step_num = retval.geometry_step.step_num
              retval.geometry_steps << GeometryStep.new 
              retval.geometry_step.step_num = last_step_num + 1
              retval.geometry_step.geometry = OutputParser.parse_updated_geometry(f, n_atoms)
            #  retval.geometry_step.geometry.lattice_vectors = vectors
              
            when /\ Final\ atomic\ structure\:/
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

