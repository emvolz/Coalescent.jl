using Coalescent
using Test
using YAML
using DataFrames
using OrdinaryDiffEq

@testset "Coalescent.jl" begin
    @testset "ModelFGY" begin
        # Test ModelFGY construction from YAML
        model_yaml = """
        modelname: TestModel
        births:
          - source: A
            recipient: B
            rate: "0.1 * A"
        deaths:
          - deme: A
            rate: "0.05 * A"
        dynamic_variables:
          - name: A
            initial_value: 100
          - name: B
            initial_value: 50
        time:
          initial: 0
          final: 10
        parameters:
          - name: beta
            value: 0.1
        """
        
        model = ModelFGY(; confstr=model_yaml)
        
        @test model.modelname == "TestModel"
        @test length(model.birthrxn) == 1
        @test length(model.deathrxn) == 1
        @test model.demes == ["A", "B"]
        @test model.initial["A"] == 100
        @test model.initial["B"] == 50
        @test model.t0 == 0
        @test model.tfin == 10
        @test model.parameters["beta"] == 0.1
    end

    @testset "solveodes" begin
        model_yaml = """
        modelname: TestModel
        births:
          - source: A
            recipient: A
            rate: "beta * A"
        deaths:
          - deme: A
            rate: "mu * A"
        dynamic_variables:
          - name: A
            initial_value: 100
        time:
          initial: 0
          final: 10
        parameters:
          - name: beta
            value: 0.1
          - name: mu
            value: 0.05
        """
        
        model = ModelFGY(; confstr=model_yaml)
        sol = solveodes(model)
        
        @test sol.t[1] == 0
        @test sol.t[end] == 10
        @test length(sol.u) > 1
        @test all(u[1] > 0 for u in sol.u)  # Population should remain positive
    end

    @testset "SampleConfiguration" begin
        sample_yaml = """
        sample:
          - time: 2
            deme: A
            size: 10
          - time: 5
            deme: B
            size: 5
        """
        
        sample_config = SampleConfiguration(; confstr=sample_yaml)
        
        @test length(sample_config.sconf) == 15
        @test count(x -> x[1] == "A" && x[2] == 0, sample_config.sconf) == 10
        @test count(x -> x[1] == "B" && x[2] == 5, sample_config.sconf) == 5
    end

    @testset "SimTree" begin
        model_yaml = """
        modelname: TestModel
        births:
          - source: A
            recipient: A
            rate: "0.1 * A"
        deaths:
          - deme: A
            rate: "0.05 * A"
        dynamic_variables:
          - name: A
            initial_value: 100
        time:
          initial: 0
          final: 10
        parameters:
          - name: beta
            value: 0.1
        """
        
        sample_yaml = """
        sample:
          - time: 10
            deme: A
            size: 10
        """
        
        model = ModelFGY(; confstr=model_yaml)
        sample_config = SampleConfiguration(; confstr=sample_yaml)
        
        tree = SimTree(model, sample_config)
        tree2 = SimTree( 10, 10 )
        
        @test tree.n == 10  # Number of samples
        @test tree2.n == 10  # Number of samples
        @test tree.nNode == 9  # n - 1 internal nodes for a binary tree
        @test tree2.nNode == 9  # n - 1 internal nodes for a binary tree
        @test tree isa SimTree 
        @test tree2 isa SimTree 
        @test length(tree.tiplabs) == 10
        @test all(!isnothing(h) for h in tree.heights)
        @test all(!isnothing(e) for e in tree.edgelength)
    end

    @testset "tonewick" begin
        model_yaml = """
        modelname: TestModel
        births:
          - source: A
            recipient: A
            rate: "0.1 * A"
        deaths:
          - deme: A
            rate: "0.05 * A"
        dynamic_variables:
          - name: A
            initial_value: 100
        time:
          initial: 0
          final: 10
        parameters:
          - name: beta
            value: 0.1
        """
        
        sample_yaml = """
        sample:
          - time: 9 
            deme: A
            size: 5
        """
        
        model = ModelFGY(; confstr=model_yaml)
        sample_config = SampleConfiguration(; confstr=sample_yaml)
        
        tree = SimTree(model, sample_config)
        newick = tonewick(tree)
        
        @test typeof(newick) == String
        @test newick[end] == ';'  # Newick strings end with a semicolon
        @test count(c -> c == '(', newick) == count(c -> c == ')', newick)  # Balanced parentheses
    end
end
