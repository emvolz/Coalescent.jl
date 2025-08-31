#= 
Test the fixed DSL implementation
=#

include("src/dsl_fixed.jl")

println("🏴‍☠️ Testing Fixed Coalescent.jl DSL")
println("=" ^ 50)

# Test 1: Basic SIR Model
println("Test 1: Basic SIR Model...")
try
    sir_model = @model "SIR_Test" begin
        @parameter β 3.0
        @parameter γ 2.0
        
        @deme I 1.0
        @nondeme S 1000.0 ode=(-β*S*I/N)
        @nondeme R 0.0 ode=(γ*I)
        
        @helper N (S + I + R)
        
        @birth (I => I) rate=(β*S*I/N)
        @death I rate=(γ*I)
        
        @timespan 1.0 35.0
    end
    
    println("✅ SIR Model created successfully!")
    println(sir_model)
    
    # Verify structure
    @assert sir_model.name == "SIR_Test"
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

# Test 2: SEIR with migrations
println("\\nTest 2: SEIR Model with Migration...")
try
    seir_model = @model "SEIR_Test" begin
        @parameter β 3.0
        @parameter γ1 2.0
        @parameter γ2 2.0
        
        @deme I 1.0
        @deme E 0.0
        @nondeme S 1e5 ode=(-β*S*I/N)
        @nondeme R 0.0 ode=(γ2*I)
        
        @helper N (S + E + I + R)
        
        @birth (I => E) rate=(β*S*I/N)
        @migration (E => I) rate=(γ1*E)
        @death I rate=(γ2*I)
        
        @timespan 1.0 35.0
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

# Test 3: Validation
println("\\nTest 3: Validation Tests...")

# Test: Model without demes should fail
try
    @model "Bad_Model" begin
        @parameter β 3.0
        @nondeme S 1000.0 ode=(-β*S)
        @timespan 0.0 10.0
    end
    println("❌ Should have failed - no demes!")
catch e
    println("✅ Correctly caught validation error: no demes")
end

# Test: Non-deme without ODE should fail
try
    @model "Bad_Model2" begin
        @parameter β 3.0
        @deme I 1.0
        @nondeme S 1000.0  # Missing ode=
        @timespan 0.0 10.0
    end
    println("❌ Should have failed - missing ODE!")
catch e
    println("✅ Correctly caught validation error: missing ODE")
end

println("\\n🌊 All DSL tests passed! The treasure chest is ready!")
println("🏴‍☠️ Yo ho ho! The DSL be seaworthy and ready for action!")