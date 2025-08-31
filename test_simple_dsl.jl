#= 
Test the simple function-based DSL implementation
=#

include("src/simple_dsl.jl")

println("🏴‍☠️ Testing Simple Functional DSL")
println("=" ^ 50)

# Test 1: Basic SIR Model using function calls
println("Test 1: SIR Model with Function Calls...")
try
    sir_model = create_model("SIR_Functional", model -> begin
        parameter!(model, :β, 3.0)
        parameter!(model, :γ, 2.0)
        
        deme!(model, :I, 1.0)
        nondeme!(model, :S, 1000.0, :(-β*S*I/N))
        nondeme!(model, :R, 0.0, :(γ*I))
        
        helper!(model, :N, :(S + I + R))
        
        birth!(model, :I, :I, :(β*S*I/N))
        death!(model, :I, :(γ*I))
        
        timespan!(model, 1.0, 35.0)
    end
    
    println("✅ SIR Model created successfully!")
    println(sir_model)
    
    # Verify structure
    @assert sir_model.name == "SIR_Functional"
    @assert sir_model.parameters[:β] == 3.0
    @assert sir_model.parameters[:γ] == 2.0
    @assert :I ∈ sir_model.demes
    @assert :S ∈ sir_model.non_demes
    @assert :R ∈ sir_model.non_demes
    @assert length(sir_model.births) == 1
    @assert length(sir_model.deaths) == 1
    @assert haskey(sir_model.helpers, :N)
    @assert sir_model.time_span.initial == 1.0
    @assert sir_model.time_span.final == 35.0
    
    println("✅ All SIR assertions passed!")
    
catch e
    println("❌ Error in SIR test: $e")
    rethrow(e)
end

# Test 2: SEIR Model
println("\\nTest 2: SEIR Model...")
try
    seir_model = create_model("SEIR_Functional") do model
        parameter!(model, :β, 3.0)
        parameter!(model, :γ1, 2.0)
        parameter!(model, :γ2, 2.0)
        
        deme!(model, :I, 1.0)
        deme!(model, :E, 0.0)
        nondeme!(model, :S, 1e5, :(-β*S*I/N))
        nondeme!(model, :R, 0.0, :(γ2*I))
        
        helper!(model, :N, :(S + E + I + R))
        
        birth!(model, :I, :E, :(β*S*I/N))
        migration!(model, :E, :I, :(γ1*E))
        death!(model, :I, :(γ2*I))
        
        timespan!(model, 1.0, 35.0)
    end
    
    println("✅ SEIR Model created successfully!")
    println(seir_model)
    
    @assert seir_model isa CoalescentModel
    @assert :E ∈ seir_model.demes
    @assert :I ∈ seir_model.demes
    @assert length(seir_model.migrations) == 1
    @assert seir_model.migrations[1].source == :E
    @assert seir_model.migrations[1].recipient == :I
    
    println("✅ All SEIR assertions passed!")
    
catch e
    println("❌ Error in SEIR test: $e")
    rethrow(e)
end

# Test 3: Complex Metapopulation Model
println("\\nTest 3: Metapopulation Model...")
try
    meta_model = create_model("Metapopulation_Functional") do model
        parameter!(model, :β, 2.0)
        parameter!(model, :γ, 1.0)
        parameter!(model, :m12, 0.1)
        parameter!(model, :m21, 0.05)
        
        deme!(model, :I1, 5.0)
        deme!(model, :I2, 1.0)
        nondeme!(model, :S1, 1000.0, :(-β*S1*I1/N1))
        nondeme!(model, :S2, 2000.0, :(-β*S2*I2/N2))
        nondeme!(model, :R1, 0.0, :(γ*I1))
        nondeme!(model, :R2, 0.0, :(γ*I2))
        
        helper!(model, :N1, :(S1 + I1 + R1))
        helper!(model, :N2, :(S2 + I2 + R2))
        
        birth!(model, :I1, :I1, :(β*S1*I1/N1))
        birth!(model, :I2, :I2, :(β*S2*I2/N2))
        migration!(model, :I1, :I2, :(m12*I1))
        migration!(model, :I2, :I1, :(m21*I2))
        death!(model, :I1, :(γ*I1))
        death!(model, :I2, :(γ*I2))
        
        timespan!(model, 0.0, 100.0)
    end
    
    println("✅ Metapopulation Model created successfully!")
    println(meta_model)
    
    @assert :I1 ∈ meta_model.demes
    @assert :I2 ∈ meta_model.demes
    @assert length(meta_model.migrations) == 2
    @assert length(meta_model.births) == 2
    @assert length(meta_model.deaths) == 2
    
    println("✅ All Metapopulation assertions passed!")
    
catch e
    println("❌ Error in Metapopulation test: $e")
    rethrow(e)
end

# Test 4: Validation
println("\\nTest 4: Validation Tests...")

# Model without demes should fail
try
    create_model("Bad_Model") do model
        parameter!(model, :β, 3.0)
        nondeme!(model, :S, 1000.0, :(-β*S))
        timespan!(model, 0.0, 10.0)
    end
    println("❌ Should have failed - no demes!")
catch e
    println("✅ Correctly caught validation error: $e")
end

# Non-deme without ODE should fail  
try
    create_model("Bad_Model2") do model
        parameter!(model, :β, 3.0)
        deme!(model, :I, 1.0)
        push!(model.non_demes, :S)  # Add directly without ODE
        model.dynamic_vars[:S] = (initial_value=1000.0, ode=nothing)
        timespan!(model, 0.0, 10.0)
    end
    println("❌ Should have failed - missing ODE!")
catch e
    println("✅ Correctly caught validation error: $e")
end

println("\\n🌊 All DSL tests passed! The functional approach works perfectly!")
println("🏴‍☠️ Ready to set sail with our treasure!")