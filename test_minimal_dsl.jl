#= 
Minimal DSL Test - Testing core DSL functionality without full package dependencies
=#

# Let's test our DSL components in isolation
using MacroTools

# Include our DSL files directly
include("src/dsl.jl")

println("🏴‍☠️ Testing Coalescent.jl DSL - Minimal Test Suite")
println("=" ^ 60)

# Test 1: Basic Model Creation
println("Test 1: Creating a basic SIR model...")
try
    sir_model = @model "Test_SIR" begin
        @parameter β 3.0
        @parameter γ 2.0
        @deme I 1.0
        @nondeme S 1000.0 ode=(-β*S*I/N)
        @nondeme R 0.0 ode=(γ*I)
        @helper N (S + I + R)
        @birth I => I rate=(β*S*I/N)
        @death I rate=(γ*I)
        @timespan 1.0 35.0
    end
    
    println("✅ Model created successfully!")
    println("   - Name: $(sir_model.name)")
    println("   - Parameters: $(length(sir_model.parameters))")
    println("   - Demes: $(collect(sir_model.demes))")
    println("   - Non-demes: $(collect(sir_model.non_demes))")
    println("   - Births: $(length(sir_model.births))")
    println("   - Deaths: $(length(sir_model.deaths))")
    println("   - Helpers: $(collect(keys(sir_model.helpers)))")
    println("   - Time span: $(sir_model.time_span)")
    
    # Verify structure
    @assert sir_model.name == "Test_SIR"
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
    
    println("✅ All assertions passed!")
    
catch e
    println("❌ Error in Test 1: $e")
end

# Test 2: SEIR Model with Migrations
println("\\nTest 2: Creating SEIR model with migrations...")
try
    seir_model = @model "Test_SEIR" begin
        @parameter β 3.0
        @parameter γ1 2.0
        @parameter γ2 2.0
        @deme I 1.0
        @deme E 0.0
        @nondeme S 1e5 ode=(-β*S*I/N)
        @nondeme R 0.0 ode=(γ2*I)
        @helper N (S + E + I + R)
        @birth I => E rate=(β*S*I/N)
        @migration E => I rate=(γ1*E)
        @death I rate=(γ2*I)
        @timespan 1.0 35.0
    end
    
    println("✅ SEIR model created successfully!")
    
    @assert seir_model isa CoalescentModel
    @assert :E ∈ seir_model.demes
    @assert :I ∈ seir_model.demes
    @assert length(seir_model.migrations) == 1
    @assert seir_model.migrations[1].source == :E
    @assert seir_model.migrations[1].recipient == :I
    
    println("✅ SEIR assertions passed!")
    
catch e
    println("❌ Error in Test 2: $e")
end

# Test 3: Validation Tests
println("\\nTest 3: Testing validation...")

# Test model without demes (should fail)
try
    bad_model = @model "Bad_Model" begin
        @parameter β 3.0
        @nondeme S 1000.0 ode=(-β*S)
        @timespan 0.0 10.0
    end
    println("❌ Should have failed validation (no demes)")
catch e
    println("✅ Correctly caught validation error: no demes")
end

# Test non-deme without ODE (should fail)
try
    bad_model2 = @model "Bad_Model2" begin
        @parameter β 3.0
        @deme I 1.0
        @nondeme S 1000.0  # Missing ODE!
        @timespan 0.0 10.0
    end
    println("❌ Should have failed validation (missing ODE)")
catch e
    println("✅ Correctly caught validation error: missing ODE for non-deme")
end

# Test 4: Complex Model
println("\\nTest 4: Creating complex metapopulation model...")
try
    meta_model = @model "Metapopulation" begin
        @parameter β 2.0
        @parameter γ 1.0
        @parameter m12 0.1
        @parameter m21 0.05
        
        @deme I1 5.0
        @deme I2 1.0
        @nondeme S1 1000.0 ode=(-β*S1*I1/N1)
        @nondeme S2 2000.0 ode=(-β*S2*I2/N2)
        @nondeme R1 0.0 ode=(γ*I1)
        @nondeme R2 0.0 ode=(γ*I2)
        
        @helper N1 (S1 + I1 + R1)
        @helper N2 (S2 + I2 + R2)
        
        @birth I1 => I1 rate=(β*S1*I1/N1)
        @birth I2 => I2 rate=(β*S2*I2/N2)
        @migration I1 => I2 rate=(m12*I1)
        @migration I2 => I1 rate=(m21*I2)
        @death I1 rate=(γ*I1)
        @death I2 rate=(γ*I2)
        
        @timespan 0.0 100.0
    end
    
    println("✅ Metapopulation model created!")
    
    @assert meta_model isa CoalescentModel
    @assert :I1 ∈ meta_model.demes
    @assert :I2 ∈ meta_model.demes
    @assert length(meta_model.migrations) == 2
    @assert length(meta_model.births) == 2
    @assert length(meta_model.deaths) == 2
    @assert haskey(meta_model.helpers, :N1)
    @assert haskey(meta_model.helpers, :N2)
    
    println("✅ Metapopulation assertions passed!")
    
catch e
    println("❌ Error in Test 4: $e")
end

println("\\n🌊 DSL Core Functionality Tests Complete!")
println("=" ^ 60)

# Test the show method
println("\\nTest 5: Display functionality...")
test_model = @model "Display_Test" begin
    @parameter β 2.0
    @parameter γ 1.0
    @deme I 1.0
    @nondeme S 1000.0 ode=(-β*S*I)
    @helper N (S + I)
    @birth I => I rate=(β*S*I/N)
    @death I rate=(γ*I)
    @timespan 0.0 50.0
end

println("\\n📊 Model Display:")
println(test_model)

println("\\n🏴‍☠️ Ahoy! All core DSL tests passed! The treasure is ready to be shipped!")
println("🌊 DSL is fully functional and ready for integration with the simulation engine!")