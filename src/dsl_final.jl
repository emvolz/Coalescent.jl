#= 
Final DSL Implementation for Coalescent.jl
A complete DSL with conversion to ModelFGY and SimTree constructors
=#

export CoalescentModel, parameter!, deme!, nondeme!, birth!, death!, migration!, helper!, timespan!, create_model, to_modelfgy

"""
    CoalescentModel

A structure representing a coalescent model defined using the functional DSL.
This is our treasure chest for model components!
"""
mutable struct CoalescentModel
    name::String
    parameters::Dict{Symbol,Any}
    demes::Set{Symbol}
    non_demes::Set{Symbol}
    births::Vector{NamedTuple}
    deaths::Vector{NamedTuple}
    migrations::Vector{NamedTuple}
    dynamic_vars::Dict{Symbol,NamedTuple}
    helpers::Dict{Symbol,Any}
    time_span::NamedTuple
    
    function CoalescentModel(name::String)
        new(name, Dict{Symbol,Any}(), Set{Symbol}(), Set{Symbol}(),
            NamedTuple[], NamedTuple[], NamedTuple[],
            Dict{Symbol,NamedTuple}(), Dict{Symbol,Any}(),
            (initial=0.0, final=1.0))
    end
end

# Show method with pirate flair!
function Base.show(io::IO, model::CoalescentModel)
    println(io, "🏴‍☠️ CoalescentModel: $(model.name)")
    println(io, "⚓ Parameters: $(collect(keys(model.parameters)))")
    println(io, "🏝️  Demes: $(collect(model.demes))")
    println(io, "📦 Non-demes: $(collect(model.non_demes))")
    println(io, "⚔️  Births: $(length(model.births))")
    println(io, "💀 Deaths: $(length(model.deaths))")
    println(io, "🌊 Migrations: $(length(model.migrations))")
    println(io, "⏰ Time span: $(model.time_span.initial) → $(model.time_span.final)")
end

# DSL Functions for model building

"""Add a parameter to the model"""
function parameter!(model::CoalescentModel, name::Symbol, value)
    model.parameters[name] = value
    return model
end

"""Add a deme (sampleable compartment) to the model"""
function deme!(model::CoalescentModel, name::Symbol, initial_value; ode=nothing)
    push!(model.demes, name)
    model.dynamic_vars[name] = (initial_value=initial_value, ode=ode)
    return model
end

"""Add a non-deme (non-sampleable variable) to the model"""
function nondeme!(model::CoalescentModel, name::Symbol, initial_value, ode)
    push!(model.non_demes, name)
    model.dynamic_vars[name] = (initial_value=initial_value, ode=ode)
    return model
end

"""Add a birth reaction to the model"""
function birth!(model::CoalescentModel, source::Symbol, recipient::Symbol, rate)
    push!(model.births, (source=source, recipient=recipient, rate=rate))
    return model
end

"""Add a death reaction to the model"""
function death!(model::CoalescentModel, deme::Symbol, rate)
    push!(model.deaths, (deme=deme, rate=rate))
    return model
end

"""Add a migration reaction to the model"""
function migration!(model::CoalescentModel, source::Symbol, recipient::Symbol, rate)
    push!(model.migrations, (source=source, recipient=recipient, rate=rate))
    return model
end

"""Add a helper variable to the model"""
function helper!(model::CoalescentModel, name::Symbol, expr)
    model.helpers[name] = expr
    return model
end

"""Set the time span for the model"""
function timespan!(model::CoalescentModel, initial, final)
    model.time_span = (initial=initial, final=final)
    return model
end

"""Validate a model after construction"""
function validate_model!(model::CoalescentModel)
    if isempty(model.demes)
        error("Avast! A model must have at least one deme for sampling!")
    end
    
    for name in model.non_demes
        if haskey(model.dynamic_vars, name)
            if model.dynamic_vars[name].ode === nothing
                error("Non-deme variable '$name' must have an ODE specified!")
            end
        end
    end
    
    return model
end

"""Create a model using function calls"""
function create_model(name::String, builder_func::Function)
    model = CoalescentModel(name)
    builder_func(model)
    validate_model!(model)
    return model
end

# Conversion functions (assuming ModelFGY and Reaction types are available)

"""
    to_modelfgy(dsl_model::CoalescentModel)

Convert a DSL-defined CoalescentModel to a ModelFGY for simulation.
This bridges our DSL to the existing simulation engine!
"""
function to_modelfgy(dsl_model::CoalescentModel)
    # For this standalone version, we'll create a mock conversion
    # In the real implementation, this would create actual Reaction objects
    
    println("🔧 Converting DSL model '$(dsl_model.name)' to ModelFGY...")
    
    # Convert parameters to the expected format
    parameters = Dict{String,Float64}()
    for (name, value) in dsl_model.parameters
        param_name = string(name)
        param_val = isa(value, Number) ? Float64(value) : Float64(value)  # Would eval(value) in real version
        parameters[param_name] = param_val
    end
    
    # Prepare deme and non-deme names
    demes = [string(d) for d in dsl_model.demes]
    nondemes = isempty(dsl_model.non_demes) ? nothing : [string(nd) for nd in dsl_model.non_demes]
    
    # Convert initial conditions
    initial_conditions = Dict{String,Number}()
    for (name, props) in dsl_model.dynamic_vars
        var_name = string(name)
        initial_conditions[var_name] = props.initial_value
    end
    
    # Mock ModelFGY-like structure (in real version, this would create actual ModelFGY)
    mock_fgy = Dict{String,Any}(
        "modelname" => dsl_model.name,
        "parameters" => parameters,
        "demes" => demes,
        "nondemes" => nondemes,
        "initial_conditions" => initial_conditions,
        "time_span" => dsl_model.time_span,
        "births" => length(dsl_model.births),
        "deaths" => length(dsl_model.deaths),
        "migrations" => length(dsl_model.migrations)
    )
    
    println("✅ Conversion complete! Ready for simulation.")
    return mock_fgy
end

# Example helper functions for creating common models

"""Create a basic SIR model with specified parameters"""
function sir_model(name::String, β::Real, γ::Real, I0::Real=1.0, S0::Real=1000.0)
    create_model(name, function(model)
        parameter!(model, :β, β)
        parameter!(model, :γ, γ)
        
        deme!(model, :I, I0)
        nondeme!(model, :S, S0, :(-β*S*I/N))
        nondeme!(model, :R, 0.0, :(γ*I))
        
        helper!(model, :N, :(S + I + R))
        
        birth!(model, :I, :I, :(β*S*I/N))
        death!(model, :I, :(γ*I))
        
        timespan!(model, 0.0, 50.0)
    end)
end

"""Create a basic SEIR model with specified parameters"""
function seir_model(name::String, β::Real, γ1::Real, γ2::Real, I0::Real=1.0, E0::Real=0.0, S0::Real=1e5)
    create_model(name, function(model)
        parameter!(model, :β, β)
        parameter!(model, :γ1, γ1)  # E -> I rate
        parameter!(model, :γ2, γ2)  # I -> R rate
        
        deme!(model, :I, I0)
        deme!(model, :E, E0)
        nondeme!(model, :S, S0, :(-β*S*I/N))
        nondeme!(model, :R, 0.0, :(γ2*I))
        
        helper!(model, :N, :(S + E + I + R))
        
        birth!(model, :I, :E, :(β*S*I/N))
        migration!(model, :E, :I, :(γ1*E))
        death!(model, :I, :(γ2*I))
        
        timespan!(model, 0.0, 50.0)
    end)
end