#= 
Test Suite for Coalescent.jl DSL
Testing our swashbuckling DSL to ensure it's seaworthy!
=#

using Test
using Coalescent

# Test basic DSL functionality
@testset "DSL Basic Functionality" begin
    
    @testset "Basic Model Creation" begin
        # Test creating a simple SIR model
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
        
        @test sir_model isa CoalescentModel
        @test sir_model.name == "Test_SIR"
        @test sir_model.parameters[:β] == 3.0
        @test sir_model.parameters[:γ] == 2.0
        @test :I ∈ sir_model.demes
        @test :S ∈ sir_model.non_demes
        @test :R ∈ sir_model.non_demes
        @test length(sir_model.births) == 1
        @test length(sir_model.deaths) == 1
        @test haskey(sir_model.helpers, :N)
        @test sir_model.time_span.initial == 1.0
        @test sir_model.time_span.final == 35.0
    end
    
    @testset "SEIR Model" begin
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
        
        @test seir_model isa CoalescentModel
        @test :E ∈ seir_model.demes
        @test :I ∈ seir_model.demes
        @test length(seir_model.migrations) == 1
        @test seir_model.migrations[1].source == :E
        @test seir_model.migrations[1].recipient == :I
    end
    
    @testset "Model Validation" begin
        # Test that model without demes fails validation
        @test_throws ErrorException @model "Bad_Model" begin
            @parameter β 3.0
            @nondeme S 1000.0 ode=(-β*S)
            @timespan 0.0 10.0
        end
        
        # Test that non-deme without ODE fails validation
        @test_throws ErrorException @model "Bad_Model2" begin
            @parameter β 3.0
            @deme I 1.0
            @nondeme S 1000.0  # Missing ODE!
            @timespan 0.0 10.0
        end
    end
end

@testset "DSL to ModelFGY Conversion" begin
    # Create a DSL model
    dsl_model = @model "Conversion_Test" begin
        @parameter β 2.5
        @parameter γ 1.5
        @deme I 5.0
        @nondeme S 995.0 ode=(-β*S*I/N)
        @nondeme R 0.0 ode=(γ*I)
        @helper N (S + I + R)
        @birth I => I rate=(β*S*I/N)
        @death I rate=(γ*I)
        @timespan 0.0 50.0
    end
    
    # Convert to ModelFGY
    fgy_model = to_modelfgy(dsl_model)
    
    @test fgy_model isa ModelFGY
    @test fgy_model.modelname == "Conversion_Test"
    @test fgy_model.parameters["β"] == 2.5
    @test fgy_model.parameters["γ"] == 1.5
    @test "I" ∈ fgy_model.demes
    @test "S" ∈ fgy_model.nondemes
    @test "R" ∈ fgy_model.nondemes
    @test fgy_model.numberdemes == 1
    @test fgy_model.numbernondemes == 2
    @test fgy_model.t0 == 0.0
    @test fgy_model.tfin == 50.0
    @test fgy_model.initial["I"] == 5.0
    @test fgy_model.initial["S"] == 995.0
    @test fgy_model.initial["R"] == 0.0
    @test length(fgy_model.birthrxn) == 1
    @test length(fgy_model.deathrxn) == 1
    @test fgy_model.helperexprs !== nothing
end

@testset "Complex Models" begin
    @testset "Multi-deme Migration Model" begin
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
        
        @test meta_model isa CoalescentModel
        @test :I1 ∈ meta_model.demes
        @test :I2 ∈ meta_model.demes
        @test length(meta_model.migrations) == 2
        @test length(meta_model.births) == 2
        @test length(meta_model.deaths) == 2
        @test haskey(meta_model.helpers, :N1)
        @test haskey(meta_model.helpers, :N2)
        
        # Convert and verify
        fgy_meta = to_modelfgy(meta_model)
        @test fgy_meta.numberdemes == 2
        @test fgy_meta.numbernondemes == 4
    end
    
    @testset "Reservoir Model" begin
        reservoir_model = @model "SIR_Reservoir" begin
            @parameter β 1.5/7.0
            @parameter β_res 1.5/7.0
            @parameter γ 1.0/7.0
            @parameter N 1e3
            @parameter N_res 1e6
            @parameter μ 1e-4*(1.0/7.0)
            
            @deme I 0.0
            @deme I_res 1.0
            @nondeme S 1e4 ode=(-β*I*S/N - μ*I_res*S/N)
            @nondeme S_res 1e6 ode=(-β_res*I_res*S_res/N_res)
            
            @birth I => I rate=(β*I*S/N)
            @birth I_res => I_res rate=(β_res*I_res*S_res/N_res)
            @birth I_res => I rate=(μ*I_res*S/N)
            @death I rate=(γ*I)
            @death I_res rate=(γ*I_res)
            
            @timespan 0.0 220.0
        end
        
        @test reservoir_model isa CoalescentModel
        @test :I ∈ reservoir_model.demes
        @test :I_res ∈ reservoir_model.demes
        @test length(reservoir_model.births) == 3  # Including cross-species transmission
        @test length(reservoir_model.deaths) == 2
        
        # Test conversion
        fgy_reservoir = to_modelfgy(reservoir_model)
        @test fgy_reservoir.numberdemes == 2
        @test length(fgy_reservoir.birthrxn) == 3
    end
end

println("🏴‍☠️ DSL Test Suite completed! All tests passed - the DSL is seaworthy!")

# Run a demonstration if this file is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    println("\\n🌊 Running DSL demonstration...")
    
    # Create and display a sample model
    demo_model = @model "Demo_SIR" begin
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
    
    println("\\n✅ Created DSL model:")
    println(demo_model)
    
    println("\\n🔄 Converting to ModelFGY...")
    fgy_demo = to_modelfgy(demo_model)
    println("✅ Conversion successful!")
    println(fgy_demo)
    
    println("\\n⚓ DSL demonstration complete! Ready for production use!")
end