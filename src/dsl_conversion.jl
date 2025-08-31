#= 
DSL to ModelFGY Conversion
Converting our DSL treasure into the format the simulation engine understands!
=#

using MacroTools

"""
    to_modelfgy(dsl_model::CoalescentModel)::ModelFGY

Convert a DSL-defined CoalescentModel to a ModelFGY for simulation.
This be the bridge between our shiny DSL and the battle-tested engine!

# Arguments
- `dsl_model::CoalescentModel`: The model defined using our DSL

# Returns
- `ModelFGY`: The converted model ready for simulation

# Example
```julia
sir_dsl = @model "SIR" begin
    # ... model definition ...
end

sir_fgy = to_modelfgy(sir_dsl)
tree = SimTree(sir_fgy, sample_config)
```
"""
function to_modelfgy(dsl_model::CoalescentModel)::ModelFGY
    # Convert parameters to the expected format
    parameters = Dict{String,Float64}()
    for (name, value) in dsl_model.parameters
        param_name = string(name)
        param_val = isa(value, Number) ? Float64(value) : Float64(eval(value))
        parameters[param_name] = param_val
    end
    
    # Convert birth reactions to Reaction objects
    birth_reactions = Reaction[]
    for birth in dsl_model.births
        push!(birth_reactions, Reaction(
            string(birth.source),
            string(birth.recipient), 
            RXN_BIRTH,
            _convert_expr_to_string_vars(birth.rate)
        ))
    end
    
    # Convert migration reactions to Reaction objects
    migration_reactions = Reaction[]
    for migration in dsl_model.migrations
        push!(migration_reactions, Reaction(
            string(migration.source),
            string(migration.recipient),
            RXN_MIG,
            _convert_expr_to_string_vars(migration.rate)
        ))
    end
    
    # Convert death reactions to Reaction objects
    death_reactions = Reaction[]
    for death in dsl_model.deaths
        push!(death_reactions, Reaction(
            string(death.deme),
            RXN_DEATH,
            _convert_expr_to_string_vars(death.rate)
        ))
    end
    
    # Prepare deme names and non-deme names
    demes = [string(d) for d in dsl_model.demes]
    nondemes = isempty(dsl_model.non_demes) ? nothing : [string(nd) for nd in dsl_model.non_demes]
    
    # Convert initial conditions
    initial_conditions = Dict{String,Number}()
    for (name, props) in dsl_model.dynamic_vars
        var_name = string(name)
        initial_conditions[var_name] = props.initial_value
    end
    
    # Convert non-deme reactions (ODEs for non-demes)
    nondeme_reactions = Reaction[]
    for name in dsl_model.non_demes
        name_str = string(name)
        if haskey(dsl_model.dynamic_vars, name)
            ode_expr = dsl_model.dynamic_vars[name].ode
            if ode_expr !== nothing
                push!(nondeme_reactions, Reaction(
                    name_str,
                    RXN_DYNVAR,
                    _convert_expr_to_string_vars(ode_expr)
                ))
            end
        end
    end
    
    # Convert helper expressions
    helper_exprs = nothing
    if !isempty(dsl_model.helpers)
        helper_exprs = Expr[]
        for (name, expr) in dsl_model.helpers
            helper_assignment = :($(Symbol(name)) = $(_convert_expr_to_string_vars(expr)))
            push!(helper_exprs, helper_assignment)
        end
    end
    
    # Create the ModelFGY
    return ModelFGY(
        dsl_model.name,                    # modelname
        birth_reactions,                   # birthrxn
        migration_reactions,               # migrationrxn
        death_reactions,                  # deathrxn
        nondeme_reactions,                # nondemerxn
        demes,                            # demes
        nondemes,                         # nondemes
        length(demes),                    # numberdemes
        length(nondemes === nothing ? String[] : nondemes), # numbernondemes
        initial_conditions,               # initial
        Float64(dsl_model.time_span.initial),  # t0
        Float64(dsl_model.time_span.final),    # tfin
        parameters,                       # parameters
        helper_exprs                      # helperexprs
    )
end

"""
    _convert_expr_to_string_vars(expr)

Internal function to convert expressions with Symbol variables to String variables.
This be the translator between Julia symbols and string-based variable names!
"""
function _convert_expr_to_string_vars(expr)
    if isa(expr, Symbol)
        return Symbol(string(expr))  # Keep as symbol but ensure it's properly formatted
    elseif isa(expr, Expr)
        # Recursively convert all symbols in the expression
        return MacroTools.postwalk(expr) do x
            if isa(x, Symbol) && !_is_julia_builtin(x)
                # Convert variable symbols to their string representation for eval
                return Symbol(string(x))
            else
                return x
            end
        end
    else
        return expr
    end
end

"""
    _is_julia_builtin(sym::Symbol)::Bool

Check if a symbol represents a Julia builtin function or constant.
We don't want to convert these to string variables!
"""
function _is_julia_builtin(sym::Symbol)::Bool
    builtin_funcs = Set([
        :+, :-, :*, :/, :^, :sqrt, :exp, :log, :sin, :cos, :tan,
        :max, :min, :abs, :clamp, :sum, :prod, :length,
        :pi, :e, :π, :ℯ, :true, :false,
        :==, :!=, :<, :>, :<=, :>=, :&, :|, :!
    ])
    return sym ∈ builtin_funcs
end

# Add new SimTree constructors that work with DSL models

"""
    SimTree(dsl_model::CoalescentModel, sample::SampleConfiguration; computedescendants = false)

Simulate a coalescent tree from a DSL-defined model and sampling configuration.
This be the new way to set sail with yer DSL models!

# Arguments
- `dsl_model::CoalescentModel`: Model defined using the DSL macros
- `sample::SampleConfiguration`: Sampling configuration

# Keywords
- `computedescendants::Bool = false`: Whether to compute descendants for each node

# Returns
- `SimTree`: A simulated coalescent tree

# Example
```julia
# Define model using DSL
sir_model = @model "SIR" begin
    @parameter β 3.0
    @parameter γ 2.0
    @deme I 1.0
    @nondeme S 1000.0 ode=(-β*S*I/N)
    @helper N (S + I)
    @birth S => I rate=(β*S*I/N)
    @death I rate=(γ*I)
    @timespan 0.0 50.0
end

# Create sampling configuration
sample_config = SampleConfiguration(confstr = \"\"\"
sample:
  - deme: I
    time: 10.0
    size: 50
\"\"\")

# Simulate tree
tree = SimTree(sir_model, sample_config)
```
"""
function SimTree(dsl_model::CoalescentModel, sample::SampleConfiguration; computedescendants = false)
    # Convert DSL model to ModelFGY
    fgy_model = to_modelfgy(dsl_model)
    
    # Use existing SimTree constructor
    return SimTree(fgy_model, sample; computedescendants = computedescendants)
end

"""
    SimTree(dsl_model::CoalescentModel, sampletimes::Array{Float64}, samplestates::Array{String}; computedescendants = false)

Simulate a coalescent tree from a DSL model with explicit sample times and states.
For when ye know exactly when and where to drop anchor!

# Arguments
- `dsl_model::CoalescentModel`: Model defined using the DSL
- `sampletimes::Array{Float64}`: Times when samples are collected
- `samplestates::Array{String}`: Demes from which samples are collected

# Keywords
- `computedescendants::Bool = false`: Whether to compute descendants for each node

# Returns
- `SimTree`: A simulated coalescent tree

# Example
```julia
sir_model = @model "SIR" begin
    # ... model definition ...
end

# Sample 20 individuals from "I" deme at time 15.0, and 10 from "I" at time 25.0
sample_times = [fill(15.0, 20); fill(25.0, 10)]
sample_states = [fill("I", 20); fill("I", 10)]

tree = SimTree(sir_model, sample_times, sample_states)
```
"""
function SimTree(dsl_model::CoalescentModel, sampletimes::Array{Float64}, samplestates::Array{String}; computedescendants = false)
    # Convert DSL model to ModelFGY
    fgy_model = to_modelfgy(dsl_model)
    
    # Use the internal _sim_markov function directly
    return _sim_markov(fgy_model, sampletimes, samplestates, computedescendants)
end