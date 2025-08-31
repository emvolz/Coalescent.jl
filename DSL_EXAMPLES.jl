#=
Complete DSL Examples for Coalescent.jl
Showing the power and elegance of our Julia DSL!

This file demonstrates how to use the functional DSL to create complex
epidemiological models that would replace YAML-based definitions.
=#

include("src/dsl_final.jl")

println("🏴‍☠️ Welcome to the Coalescent.jl DSL Examples!")
println("Ahoy! Prepare to witness the treasure of pure Julia model definition!")
println("=" ^ 70)

# Example 1: Basic SIR Model
println("\\n⚔️  Example 1: Basic SIR Model")
println("-" ^ 40)

sir_basic = create_model("Basic_SIR", function(model)
    # Parameters - the constants that guide our voyage
    parameter!(model, :β, 3.0)      # Transmission rate
    parameter!(model, :γ, 2.0)      # Recovery rate
    
    # Demes (samplingable compartments)
    deme!(model, :I, 1.0)           # Infected individuals
    
    # Non-demes (important variables but not samplingable) 
    nondeme!(model, :S, 1000.0, :(-β*S*I/N))  # Susceptible population
    nondeme!(model, :R, 0.0, :(γ*I))          # Recovered population
    
    # Helper variables
    helper!(model, :N, :(S + I + R))          # Total population
    
    # Reactions
    birth!(model, :I, :I, :(β*S*I/N))         # Transmission
    death!(model, :I, :(γ*I))                 # Recovery
    
    # Time span
    timespan!(model, 1.0, 35.0)
end)

println(sir_basic)
println("✅ Basic SIR model created - classic epidemiology!")

# Example 2: SEIR Model with Incubation Period
println("\\n🦠 Example 2: SEIR Model (with incubation period)")
println("-" ^ 50)

seir_advanced = create_model("SEIR_Advanced", function(model)
    parameter!(model, :β, 3.0)      # Transmission rate
    parameter!(model, :σ, 2.0)      # Incubation rate (E->I)
    parameter!(model, :γ, 1.5)      # Recovery rate (I->R)
    
    # Multiple demes for sampling
    deme!(model, :I, 1.0)           # Infected (symptomatic)
    deme!(model, :E, 0.0)           # Exposed (incubating)
    
    # Non-deme compartments
    nondeme!(model, :S, 1e5, :(-β*S*I/N))     # Susceptible
    nondeme!(model, :R, 0.0, :(γ*I))          # Recovered
    
    helper!(model, :N, :(S + E + I + R))
    
    # Disease progression chain: S -> E -> I -> R
    birth!(model, :I, :E, :(β*S*I/N))         # Transmission
    migration!(model, :E, :I, :(σ*E))         # Disease progression
    death!(model, :I, :(γ*I))                 # Recovery
    
    timespan!(model, 1.0, 50.0)
end)

println(seir_advanced)
println("✅ SEIR model created - more realistic with incubation!")

# Example 3: Multi-host Model with Reservoir
println("\\n🐾 Example 3: SIR with Animal Reservoir")
println("-" ^ 45)

reservoir_model = create_model("SIR_Reservoir", function(model)
    # Human population parameters
    parameter!(model, :β_h, 1.5/7.0)         # Human transmission rate
    parameter!(model, :γ, 1.0/7.0)           # Recovery rate (both species)
    parameter!(model, :N_h, 1e3)             # Human population size
    
    # Reservoir population parameters
    parameter!(model, :β_r, 1.5/7.0)         # Reservoir transmission rate
    parameter!(model, :N_r, 1e6)             # Reservoir population size
    
    # Cross-species transmission
    parameter!(model, :μ, 1e-4*(1.0/7.0))    # Spillover rate from reservoir
    
    # Demes for both human and reservoir populations
    deme!(model, :I_h, 0.0)                  # Infected humans
    deme!(model, :I_r, 1.0)                  # Infected reservoir animals
    
    # Susceptible populations (non-samplingable)
    nondeme!(model, :S_h, 1e4, :(-β_h*I_h*S_h/N_h - μ*I_r*S_h/N_h))
    nondeme!(model, :S_r, 1e6, :(-β_r*I_r*S_r/N_r))
    
    # Within-species transmission
    birth!(model, :I_h, :I_h, :(β_h*I_h*S_h/N_h))        # Human-to-human
    birth!(model, :I_r, :I_r, :(β_r*I_r*S_r/N_r))        # Reservoir circulation
    
    # Cross-species spillover (the dangerous part!)
    birth!(model, :I_r, :I_h, :(μ*I_r*S_h/N_h))          # Reservoir to human
    
    # Recovery in both populations
    death!(model, :I_h, :(γ*I_h))                         # Human recovery
    death!(model, :I_r, :(γ*I_r))                         # Reservoir recovery
    
    timespan!(model, 0.0, 220.0)  # Longer time for reservoir dynamics
end)

println(reservoir_model)
println("✅ Reservoir model created - modeling zoonotic spillover!")

# Example 4: Metapopulation Model (Spatial Structure)
println("\\n🗺️  Example 4: Metapopulation Model")
println("-" ^ 40)

metapop_model = create_model("Metapopulation_SIR", function(model)
    # Disease parameters (same for both patches)
    parameter!(model, :β, 2.5)
    parameter!(model, :γ, 1.0)
    
    # Migration parameters between patches
    parameter!(model, :m12, 0.1)      # Migration rate patch 1 -> 2
    parameter!(model, :m21, 0.05)     # Migration rate patch 2 -> 1
    
    # Patch 1 (urban center)
    deme!(model, :I1, 5.0)            # Initial outbreak
    nondeme!(model, :S1, 1000.0, :(-β*S1*I1/N1))
    nondeme!(model, :R1, 0.0, :(γ*I1))
    
    # Patch 2 (rural area)  
    deme!(model, :I2, 1.0)            # Smaller initial cases
    nondeme!(model, :S2, 2000.0, :(-β*S2*I2/N2))
    nondeme!(model, :R2, 0.0, :(γ*I2))
    
    # Patch-specific population sizes
    helper!(model, :N1, :(S1 + I1 + R1))
    helper!(model, :N2, :(S2 + I2 + R2))
    
    # Within-patch transmission
    birth!(model, :I1, :I1, :(β*S1*I1/N1))
    birth!(model, :I2, :I2, :(β*S2*I2/N2))
    
    # Between-patch migration (infected individuals moving)
    migration!(model, :I1, :I2, :(m12*I1))     # Urban to rural
    migration!(model, :I2, :I1, :(m21*I2))     # Rural to urban
    
    # Recovery in both patches
    death!(model, :I1, :(γ*I1))
    death!(model, :I2, :(γ*I2))
    
    timespan!(model, 0.0, 100.0)
end)

println(metapop_model)
println("✅ Metapopulation model created - spatial epidemiology!")

# Example 5: Using Helper Functions for Quick Model Creation
println("\\n🚀 Example 5: Quick Model Creation with Helpers")
println("-" ^ 50)

# Create models using the helper functions
quick_sir = sir_model("Quick_SIR", 2.5, 1.0, 5.0, 2000.0)
quick_seir = seir_model("Quick_SEIR", 3.0, 2.0, 1.5, 1.0, 0.0, 5000.0)

println("SIR Helper Result:")
println(quick_sir)
println("\\nSEIR Helper Result:")  
println(quick_seir)
println("✅ Helper functions make model creation blazingly fast!")

# Example 6: Model Comparison and Conversion
println("\\n🔄 Example 6: Model Conversion to ModelFGY")
println("-" ^ 50)

# Convert our models to show how they bridge to the simulation engine
println("Converting SIR model...")
sir_fgy = to_modelfgy(sir_basic)
println("FGY Parameters: $(sir_fgy["parameters"])")
println("FGY Demes: $(sir_fgy["demes"])")
println("FGY Initial Conditions: $(sir_fgy["initial_conditions"])")

println("\\n🌊 All examples completed successfully!")
println("=" ^ 70)
println("🏴‍☠️ The DSL treasure is ready for production use!")
println("⚓ Key advantages over YAML:")
println("   • Type safety and error checking at definition time")
println("   • Native Julia expressions and functions")
println("   • Better IDE support and debugging")
println("   • Programmatic model construction and modification")
println("   • More concise and readable syntax")
println("   • Full integration with Julia's ecosystem")
println("\\n🚢 Set sail and simulate with confidence, ye savvy coder!")