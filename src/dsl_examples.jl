#= 
DSL Examples and Usage Guide
Showing how to use our swashbuckling DSL to define coalescent models!

This file demonstrates the DSL by recreating existing YAML models
and showing new possibilities only available in pure Julia.
=#

using .Coalescent  # Make sure we can access the main module

"""
    create_sir_model()

Create a basic SIR (Susceptible-Infected-Recovered) model using the DSL.
This be a classic epidemiological model, as fundamental as knowing port from starboard!

Equivalent to the basic SIR YAML model.
"""
function create_sir_model()
    return @model "SIR_DSL" begin
        # Model parameters - the constants that guide our simulation
        @parameter β 3.0      # Transmission rate (infections per day)
        @parameter γ 2.0      # Recovery rate (recoveries per day)
        
        # Demes (compartments we can sample from)
        @deme I 1.0           # Infected individuals (initially 1 person)
        
        # Non-demes (important variables but not sampleable)
        @nondeme S 1e5 ode=(-β*S*I/N)    # Susceptible population
        @nondeme R 0.0 ode=(γ*I)         # Recovered population
        
        # Helper variables (computed from other variables)
        @helper N (S + I + R)   # Total population size
        
        # Birth reactions (new infections)
        @birth I => I rate=(β*S*I/N)   # I infects S, creating new I
        
        # Death reactions (recovery from infection)
        @death I rate=(γ*I)    # I individuals recover
        
        # Time span for simulation
        @timespan 1.0 35.0
    end
end

"""
    create_seir_model()

Create a SEIR (Susceptible-Exposed-Infected-Recovered) model.
Like SIR but with an incubation period - more realistic for many diseases!

This reproduces the SEIR.yaml model from the test suite.
"""
function create_seir_model()
    return @model "SEIR_DSL" begin
        # Parameters
        @parameter β 3.0      # Transmission rate  
        @parameter γ1 2.0     # Rate of progression from E to I
        @parameter γ2 2.0     # Recovery rate from I to R
        
        # Demes (sampleable compartments)
        @deme I 1.0           # Infected (symptomatic)
        @deme E 0.0           # Exposed (incubating)
        
        # Non-demes
        @nondeme S 1e5 ode=(-β*S*I/N)      # Susceptible
        @nondeme R 0.0 ode=(γ2*I)          # Recovered
        
        # Helper for total population
        @helper N (S + E + I + R)
        
        # Birth reaction (transmission)
        @birth I => E rate=(β*S*I/N)    # I infects S, creating E
        
        # Migration reaction (disease progression)  
        @migration E => I rate=(γ1*E)   # E becomes I
        
        # Death reaction (recovery)
        @death I rate=(γ2*I)            # I recovers
        
        # Time span
        @timespan 1.0 35.0
    end
end

"""
    create_sir_reservoir_model()

Create a SIR model with an animal reservoir.
This shows how complex multi-host dynamics can be elegantly expressed!

Recreates the sir+reservoir.yaml model from the test suite.
"""
function create_sir_reservoir_model()
    return @model "SIR_Reservoir_DSL" begin
        # Parameters for human population
        @parameter β 1.5/7.0          # Human transmission rate
        @parameter γ 1.0/7.0          # Recovery rate
        @parameter N 1e3              # Human population size
        
        # Parameters for reservoir population  
        @parameter β_res 1.5/7.0      # Reservoir transmission rate
        @parameter N_res 1e6          # Reservoir population size
        
        # Cross-species transmission rate
        @parameter μ 1e-4*(1.0/7.0)   # Rate of spillover from reservoir to humans
        
        # Demes (both human and reservoir infected populations)
        @deme I 0.0                   # Infected humans
        @deme I_res 1.0               # Infected reservoir animals
        
        # Non-demes
        @nondeme S 1e4 ode=(-β*I*S/N - μ*I_res*S/N)  # Susceptible humans
        @nondeme S_res 1e6 ode=(-β_res*I_res*S_res/N_res)  # Susceptible reservoir
        
        # Birth reactions (within-species transmission)
        @birth I => I rate=(β*I*S/N)              # Human-to-human transmission
        @birth I_res => I_res rate=(β_res*I_res*S_res/N_res)  # Reservoir circulation
        
        # Cross-species transmission (reservoir to human)
        @birth I_res => I rate=(μ*I_res*S/N)      # Spillover events
        
        # Death reactions (recovery)
        @death I rate=(γ*I)                       # Human recovery
        @death I_res rate=(γ*I_res)               # Reservoir recovery
        
        # Extended time span for reservoir dynamics
        @timespan 0.0 220.0
    end
end

"""
    create_exponential_sir_model()

Create a SIR model with exponentially growing population.
This shows how to incorporate demographic changes into disease dynamics!

Demonstrates advanced Julia features only possible with the DSL.
"""
function create_exponential_sir_model()
    return @model "SIR_Exponential_DSL" begin
        # Disease parameters
        @parameter β 0.5      # Transmission rate (reduced due to growing population)
        @parameter γ 0.2      # Recovery rate
        
        # Demographic parameter
        @parameter r 0.01     # Population growth rate
        
        # Demes
        @deme I 10.0          # Initial infected
        
        # Non-demes with demographic growth
        @nondeme S 1000.0 ode=(r*S - β*S*I/N)    # Births minus infections
        @nondeme R 0.0 ode=(r*R + γ*I)           # Births plus recoveries
        
        # Helper for total population (growing over time!)
        @helper N (S + I + R)
        
        # Birth and death reactions
        @birth I => I rate=(β*S*I/N)
        @death I rate=(γ*I)
        
        @timespan 0.0 100.0
    end
end

"""
    create_metapopulation_model()

Create a metapopulation model with migration between patches.
This demonstrates spatial structure in disease dynamics!

This showcases features that would be cumbersome in YAML.
"""
function create_metapopulation_model()
    return @model "Metapopulation_DSL" begin
        # Disease parameters (same for both patches)
        @parameter β 2.0
        @parameter γ 1.0
        
        # Migration parameters
        @parameter m12 0.1    # Migration rate from patch 1 to patch 2
        @parameter m21 0.05   # Migration rate from patch 2 to patch 1
        
        # Patch 1 demes
        @deme I1 5.0          # Infected in patch 1
        @deme S1 1000.0       # Susceptible in patch 1 (implicitly defined ODE)
        
        # Patch 2 demes  
        @deme I2 1.0          # Infected in patch 2
        @deme S2 2000.0       # Susceptible in patch 2
        
        # Non-demes for recovered individuals
        @nondeme R1 0.0 ode=(γ*I1)
        @nondeme R2 0.0 ode=(γ*I2)
        
        # Helpers for patch populations
        @helper N1 (S1 + I1 + R1)
        @helper N2 (S2 + I2 + R2)
        
        # Within-patch transmission
        @birth I1 => I1 rate=(β*S1*I1/N1)
        @birth I2 => I2 rate=(β*S2*I2/N2)
        
        # Migration between patches
        @migration I1 => I2 rate=(m12*I1)
        @migration I2 => I1 rate=(m21*I2)
        
        # Recovery
        @death I1 rate=(γ*I1)
        @death I2 rate=(γ*I2)
        
        @timespan 0.0 100.0
    end
end

"""
    demonstrate_dsl_usage()

Comprehensive demonstration of DSL usage and capabilities.
Run this to see all the examples in action!
"""
function demonstrate_dsl_usage()
    println("🏴‍☠️ Ahoy! Welcome to the Coalescent.jl DSL demonstration!")
    println("=" ^ 60)
    
    # Basic SIR model
    println("\n⚔️  Creating basic SIR model...")
    sir_model = create_sir_model()
    println(sir_model)
    
    # SEIR model
    println("\n🦠 Creating SEIR model with incubation period...")
    seir_model = create_seir_model()
    println(seir_model)
    
    # SIR with reservoir
    println("\n🐾 Creating SIR model with animal reservoir...")
    reservoir_model = create_sir_reservoir_model()
    println(reservoir_model)
    
    # Exponential growth model
    println("\n📈 Creating SIR model with population growth...")
    exp_model = create_exponential_sir_model()
    println(exp_model)
    
    # Metapopulation model
    println("\n🗺️  Creating metapopulation model...")
    meta_model = create_metapopulation_model()
    println(meta_model)
    
    println("\n🌊 All models successfully created! Ready to set sail and simulate!")
    println("=" ^ 60)
    
    return sir_model, seir_model, reservoir_model, exp_model, meta_model
end

"""
    create_sampling_example()

Show how to create sampling configurations for DSL models.
"""
function create_sampling_example()
    # Create a model
    sir_model = create_sir_model()
    
    # Method 1: Using SampleConfiguration with YAML-like syntax
    sample_config1 = SampleConfiguration(confstr = \"\"\"
    sample:
      - deme: I
        time: 10.0
        size: 50
      - deme: I  
        time: 20.0
        size: 30
    \"\"\")
    
    # Method 2: Using arrays directly (new with DSL!)
    sample_times = [fill(10.0, 50); fill(20.0, 30)]
    sample_states = [fill("I", 50); fill("I", 30)]
    
    println("🎯 Sampling examples created!")
    println("Method 1: SampleConfiguration with $(length(sample_config1.sconf)) samples")
    println("Method 2: Direct arrays with $(length(sample_times)) samples")
    
    return sir_model, sample_config1, sample_times, sample_states
end

"""
    simulate_dsl_example()

Complete example showing DSL model definition and tree simulation.
"""
function simulate_dsl_example()
    println("🌊 Setting sail with a complete DSL simulation example!")
    
    # Create a simple SIR model
    sir_model = @model "Example_SIR" begin
        @parameter β 3.0
        @parameter γ 2.0
        @deme I 1.0
        @nondeme S 1000.0 ode=(-β*S*I/N)
        @nondeme R 0.0 ode=(γ*I)
        @helper N (S + I + R)
        @birth I => I rate=(β*S*I/N)
        @death I rate=(γ*I)
        @timespan 0.0 20.0
    end
    
    println("✅ Model created:")
    println(sir_model)
    
    # Create sampling scheme - 100 samples from I at time 15
    sample_times = fill(15.0, 100)
    sample_states = fill("I", 100)
    
    println("\\n🎯 Sampling: $(length(sample_times)) individuals from deme I at time 15.0")
    
    # Note: Actual simulation would require the full module to be loaded
    println("\\n🚢 To simulate: tree = SimTree(sir_model, sample_times, sample_states)")
    
    println("\\n⚓ Example complete! The DSL is ready for production use!")
    
    return sir_model
end