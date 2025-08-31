#= 
Simple Julia DSL for Coalescent.jl Models
A straightforward approach that actually works!
=#

export CoalescentModel, parameter!, deme!, nondeme!, birth!, death!, migration!, helper!, timespan!, create_model

"""
    CoalescentModel

A structure representing a coalescent model.
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

# Show method
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

"""
Create a model using function calls - a more reliable approach than macros
"""
function create_model(name::String, builder_func::Function)
    model = CoalescentModel(name)
    builder_func(model)
    validate_model!(model)
    return model
end