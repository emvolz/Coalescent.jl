#= 
Julia DSL for Coalescent.jl Models
A swashbuckling approach to defining coalescent simulation models!

This DSL provides a more idiomatic Julia approach to model definition,
replacing the YAML-based workflow with macros and native Julia syntax.
=#

using MacroTools

# Export the main DSL components
export @model, @dsl_model, CoalescentModel, SimTree

"""
    CoalescentModel

A structure representing a coalescent model defined via the DSL.
This is the treasure chest that holds all our model components!

# Fields
- `name::String`: Name of the model (like naming yer ship!)
- `parameters::Dict{Symbol,Any}`: Model parameters (the gold in our hold)
- `demes::Set{Symbol}`: Dynamic variables that can be sampled (our crew members)
- `non_demes::Set{Symbol}`: Non-deme dynamic variables (the cargo)
- `births::Vector{NamedTuple}`: Birth reactions (new recruits joining)
- `deaths::Vector{NamedTuple}`: Death reactions (sending souls to Davy Jones)
- `migrations::Vector{NamedTuple}`: Migration reactions (crew transfers)
- `dynamic_vars::Dict{Symbol,NamedTuple}`: All dynamic variables with initial conditions
- `helpers::Dict{Symbol,Expr}`: Helper expressions (navigation tools)
- `time_span::NamedTuple`: Time range (duration of our voyage)
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
    helpers::Dict{Symbol,Expr}
    time_span::NamedTuple
    
    # Constructor to initialize an empty model (like launching a new ship)
    function CoalescentModel(name::String)
        new(name, Dict{Symbol,Any}(), Set{Symbol}(), Set{Symbol}(),
            NamedTuple[], NamedTuple[], NamedTuple[],
            Dict{Symbol,NamedTuple}(), Dict{Symbol,Expr}(),
            (initial=0.0, final=1.0))
    end
end

# Show method for CoalescentModel (display our ship's manifest)
function Base.show(io::IO, model::CoalescentModel)
    println(io, "🏴‍☠️ CoalescentModel: $(model.name)")
    println(io, "⚓ Parameters: $(length(model.parameters))")
    println(io, "🏝️  Demes: $(collect(model.demes))")
    println(io, "📦 Non-demes: $(collect(model.non_demes))")
    println(io, "⚔️  Births: $(length(model.births))")
    println(io, "💀 Deaths: $(length(model.deaths))")
    println(io, "🌊 Migrations: $(length(model.migrations))")
    println(io, "⏰ Time span: $(model.time_span.initial) → $(model.time_span.final)")
end

"""
    @parameter(name, value)

Define a parameter for the model. These are the constants that guide our voyage!

# Arguments
- `name`: Parameter name (as Symbol)
- `value`: Parameter value (Number or expression)

# Example
```julia
@parameter β 2.0
@parameter γ 1.0/7.0
```
"""
macro parameter(name, value)
    return quote
        if !@isdefined(___coalescent_model___)
            error("@parameter must be used inside a @model block, ye landlubber!")
        end
        ___coalescent_model___.parameters[$(QuoteNode(name))] = $(esc(value))
    end
end

"""
    @deme(name, initial_value, [ode=nothing])

Define a deme (samplingable population compartment). 
These are the crew members we can sample from!

# Arguments
- `name`: Deme name (as Symbol)
- `initial_value`: Initial population size
- `ode`: Optional ODE for time evolution (if not specified, derived from reactions)

# Example
```julia
@deme I 1.0
@deme S 1000.0 ode=(-β*S*I/N)
```
"""
macro deme(name, initial_value, ode_expr=nothing)
    return quote
        if !@isdefined(___coalescent_model___)
            error("@deme must be used inside a @model block, savvy?")
        end
        push!(___coalescent_model___.demes, $(QuoteNode(name)))
        ___coalescent_model___.dynamic_vars[$(QuoteNode(name))] = (
            initial_value = $(esc(initial_value)),
            ode = $(ode_expr == nothing ? nothing : esc(ode_expr))
        )
    end
end

"""
    @nondeme(name, initial_value, ode)

Define a non-deme dynamic variable. These are like cargo - important but not samplingable!

# Arguments
- `name`: Variable name (as Symbol)
- `initial_value`: Initial value
- `ode`: ODE expression for time evolution (required for non-demes)

# Example
```julia
@nondeme R 0.0 ode=(γ*I)
```
"""
macro nondeme(name, initial_value, ode_expr=nothing)
    # Handle the case where ode= is specified
    if ode_expr === nothing
        return quote
            if !@isdefined(___coalescent_model___)
                error("@nondeme must be used inside a @model block, ye scallywag!")
            end
            error("Non-deme variable '$($(QuoteNode(name)))' must have an ODE specified, ye landlubber! Use: @nondeme $($(QuoteNode(name))) $($(esc(initial_value))) ode=(expression)")
        end
    else
        return quote
            if !@isdefined(___coalescent_model___)
                error("@nondeme must be used inside a @model block, ye scallywag!")
            end
            push!(___coalescent_model___.non_demes, $(QuoteNode(name)))
            ___coalescent_model___.dynamic_vars[$(QuoteNode(name))] = (
                initial_value = $(esc(initial_value)),
                ode = $(esc(ode_expr))
            )
        end
    end
end

"""
    @birth(source => recipient, rate)

Define a birth reaction. New souls joining the crew from existing members!

# Arguments
- `source => recipient`: Source and recipient demes
- `rate`: Birth rate expression

# Example
```julia
@birth I => I rate=(β*S*I/N)
@birth I => E rate=(β*S*I/N)  # Transmission from I to E
```
"""
macro birth(transfer_expr, rate_expr)
    # Parse the source => recipient expression
    if @capture(transfer_expr, source_ => recipient_)
        return quote
            if !@isdefined(___coalescent_model___)
                error("@birth must be used inside a @model block, matey!")
            end
            push!(___coalescent_model___.births, (
                source = $(QuoteNode(source)),
                recipient = $(QuoteNode(recipient)),
                rate = $(esc(rate_expr))
            ))
        end
    else
        error("@birth requires format: source => recipient, rate=(expression)")
    end
end

"""
    @death(deme, rate)

Define a death reaction. Sending souls to Davy Jones' locker!

# Arguments
- `deme`: Deme where deaths occur
- `rate`: Death rate expression

# Example
```julia
@death I rate=(γ*I)
```
"""
macro death(deme, rate_expr)
    return quote
        if !@isdefined(___coalescent_model___)
            error("@death must be used inside a @model block, ye scurvy dog!")
        end
        push!(___coalescent_model___.deaths, (
            deme = $(QuoteNode(deme)),
            rate = $(esc(rate_expr))
        ))
    end
end

"""
    @migration(source => recipient, rate)

Define a migration reaction. Crew members transferring between ships!

# Arguments
- `source => recipient`: Source and recipient demes
- `rate`: Migration rate expression

# Example
```julia
@migration E => I rate=(σ*E)  # E compartment moving to I
```
"""
macro migration(transfer_expr, rate_expr)
    # Parse the source => recipient expression
    if @capture(transfer_expr, source_ => recipient_)
        return quote
            if !@isdefined(___coalescent_model___)
                error("@migration must be used inside a @model block, ye barnacle!")
            end
            push!(___coalescent_model___.migrations, (
                source = $(QuoteNode(source)),
                recipient = $(QuoteNode(recipient)),
                rate = $(esc(rate_expr))
            ))
        end
    else
        error("@migration requires format: source => recipient, rate=(expression)")
    end
end

"""
    @helper(name, expression)

Define a helper variable. These are like navigation tools - derived from other variables!

# Arguments
- `name`: Helper variable name
- `expression`: Expression to compute the helper

# Example
```julia
@helper N (S + I + R)
```
"""
macro helper(name, expr)
    return quote
        if !@isdefined(___coalescent_model___)
            error("@helper must be used inside a @model block, ye sea dog!")
        end
        ___coalescent_model___.helpers[$(QuoteNode(name))] = $(esc(expr))
    end
end

"""
    @timespan(initial, final)

Set the time span for the simulation. The duration of our voyage!

# Arguments
- `initial`: Start time
- `final`: End time

# Example
```julia
@timespan 0.0 100.0
```
"""
macro timespan(initial, final)
    return quote
        if !@isdefined(___coalescent_model___)
            error("@timespan must be used inside a @model block, ye landlubber!")
        end
        ___coalescent_model___.time_span = (initial = $(esc(initial)), final = $(esc(final)))
    end
end

"""
    @model(name, body)

The main macro for defining a coalescent model. This is where ye chart yer course!

# Arguments
- `name`: Model name (String)
- `body`: Block containing model definition using DSL macros

# Example
```julia
sir_model = @model "SIR" begin
    # Parameters (the gold in our treasure chest)
    @parameter β 3.0
    @parameter γ 2.0
    
    # Demes (our crew members we can sample)
    @deme I 1.0
    @deme S 1e5 ode=(-β*S*I/N)
    
    # Non-demes (important cargo but not samplingable)
    @nondeme R 0.0 ode=(γ*I)
    
    # Helper variables (navigation tools)
    @helper N (S + I + R)
    
    # Birth reactions (new recruits)
    @birth S => I rate=(β*S*I/N)
    
    # Death reactions (to Davy Jones' locker)
    @death I rate=(γ*I)
    
    # Time span (duration of voyage)
    @timespan 0.0 50.0
end
```
"""
macro model(name, body)
    return quote
        # Create the model instance in our local scope
        ___coalescent_model___ = CoalescentModel($(esc(name)))
        
        # Execute the model definition body
        $(esc(body))
        
        # Validate and finalize the model
        _validate_model!(___coalescent_model___)
        
        # Return the completed model
        ___coalescent_model___
    end
end

"""
    _validate_model!(model::CoalescentModel)

Internal function to validate a model after construction.
Ensures our ship is seaworthy before setting sail!
"""
function _validate_model!(model::CoalescentModel)
    # Check that we have at least one deme for sampling
    if isempty(model.demes)
        error("Avast! A model must have at least one deme for sampling, ye scurvy dog!")
    end
    
    # Check that all referenced variables in expressions exist
    all_vars = union(model.demes, model.non_demes, keys(model.parameters), keys(model.helpers))
    
    # Validate birth reactions
    for birth in model.births
        if birth.source ∉ model.demes
            @warn "Birth source '$(birth.source)' not defined as a deme. Batten down the hatches!"
        end
        if birth.recipient ∉ model.demes
            @warn "Birth recipient '$(birth.recipient)' not defined as a deme. Check yer compass!"
        end
    end
    
    # Validate death reactions
    for death in model.deaths
        if death.deme ∉ model.demes
            @warn "Death deme '$(death.deme)' not defined. Something's fishy here!"
        end
    end
    
    # Validate migration reactions
    for migration in model.migrations
        if migration.source ∉ model.demes
            @warn "Migration source '$(migration.source)' not defined as a deme. Lost at sea!"
        end
        if migration.recipient ∉ model.demes
            @warn "Migration recipient '$(migration.recipient)' not defined as a deme. Navigate better!"
        end
    end
    
    # Ensure non-demes have ODEs specified
    for (name, props) in model.dynamic_vars
        if name ∈ model.non_demes && props.ode === nothing
            error("Non-deme variable '$name' must have an ODE specified, ye landlubber!")
        end
    end
    
    return model
end