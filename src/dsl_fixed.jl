#= 
Julia DSL for Coalescent.jl Models (Fixed Version)
A swashbuckling approach to defining coalescent simulation models!
=#

using MacroTools

# Export the main DSL components
export @model, CoalescentModel

"""
    CoalescentModel

A structure representing a coalescent model defined via the DSL.
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
    
    function CoalescentModel(name::String)
        new(name, Dict{Symbol,Any}(), Set{Symbol}(), Set{Symbol}(),
            NamedTuple[], NamedTuple[], NamedTuple[],
            Dict{Symbol,NamedTuple}(), Dict{Symbol,Expr}(),
            (initial=0.0, final=1.0))
    end
end

# Show method
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

# Helper function to parse keyword-style macro calls
function _parse_macro_call(expr)
    if @capture(expr, f_(args__))
        return f, args
    elseif @capture(expr, f_(args__, kw_args__; kwargs__))
        return f, vcat(args, kw_args, kwargs)
    else
        return nothing, [expr]
    end
end

"""
Main model macro that creates and manages the CoalescentModel
"""
macro model(name, body)
    # Parse the body and extract all the macro calls
    model_code = quote
        ___coalescent_model___ = CoalescentModel($(esc(name)))
        
        # Process each expression in the body
        $(esc(_process_model_body(body)))
        
        # Validate and return the model
        _validate_model!(___coalescent_model___)
        ___coalescent_model___
    end
    
    return model_code
end

function _process_model_body(body)
    if body.head == :block
        processed_exprs = []
        
        for expr in body.args
            if isa(expr, LineNumberNode)
                continue  # Skip line number nodes
            end
            
            processed_expr = _process_single_expr(expr)
            if processed_expr !== nothing
                push!(processed_exprs, processed_expr)
            end
        end
        
        return Expr(:block, processed_exprs...)
    else
        return _process_single_expr(body)
    end
end

function _process_single_expr(expr)
    if !isa(expr, Expr)
        return nothing
    end
    
    if @capture(expr, @parameter name_ value_)
        return quote
            ___coalescent_model___.parameters[$(QuoteNode(name))] = $(esc(value))
        end
        
    elseif @capture(expr, @deme name_ initial_value_)
        return quote
            push!(___coalescent_model___.demes, $(QuoteNode(name)))
            ___coalescent_model___.dynamic_vars[$(QuoteNode(name))] = (
                initial_value = $(esc(initial_value)),
                ode = nothing
            )
        end
        
    elseif @capture(expr, @deme name_ initial_value_ ode = ode_expr_)
        return quote
            push!(___coalescent_model___.demes, $(QuoteNode(name)))
            ___coalescent_model___.dynamic_vars[$(QuoteNode(name))] = (
                initial_value = $(esc(initial_value)),
                ode = $(esc(ode_expr))
            )
        end
        
    elseif @capture(expr, @nondeme name_ initial_value_ ode = ode_expr_)
        return quote
            push!(___coalescent_model___.non_demes, $(QuoteNode(name)))
            ___coalescent_model___.dynamic_vars[$(QuoteNode(name))] = (
                initial_value = $(esc(initial_value)),
                ode = $(esc(ode_expr))
            )
        end
        
    elseif @capture(expr, @birth (source_ => recipient_) rate = rate_expr_)
        return quote
            push!(___coalescent_model___.births, (
                source = $(QuoteNode(source)),
                recipient = $(QuoteNode(recipient)),
                rate = $(esc(rate_expr))
            ))
        end
        
    elseif @capture(expr, @death deme_ rate = rate_expr_)
        return quote
            push!(___coalescent_model___.deaths, (
                deme = $(QuoteNode(deme)),
                rate = $(esc(rate_expr))
            ))
        end
        
    elseif @capture(expr, @migration (source_ => recipient_) rate = rate_expr_)
        return quote
            push!(___coalescent_model___.migrations, (
                source = $(QuoteNode(source)),
                recipient = $(QuoteNode(recipient)),
                rate = $(esc(rate_expr))
            ))
        end
        
    elseif @capture(expr, @helper name_ helper_expr_)
        return quote
            ___coalescent_model___.helpers[$(QuoteNode(name))] = $(esc(helper_expr))
        end
        
    elseif @capture(expr, @timespan initial_ final_)
        return quote
            ___coalescent_model___.time_span = (initial = $(esc(initial)), final = $(esc(final)))
        end
        
    else
        @warn "Unrecognized expression in model body: $expr"
        return nothing
    end
end

"""
Validate a model after construction
"""
function _validate_model!(model::CoalescentModel)
    # Check that we have at least one deme
    if isempty(model.demes)
        error("Avast! A model must have at least one deme for sampling!")
    end
    
    # Check that non-demes have ODEs
    for name in model.non_demes
        if haskey(model.dynamic_vars, name)
            if model.dynamic_vars[name].ode === nothing
                error("Non-deme variable '$name' must have an ODE specified!")
            end
        end
    end
    
    return model
end