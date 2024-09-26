#= 
Structured FGY model 
=#

using YAML 
using DataFrames
using OrdinaryDiffEq
using Plots

const RXN_BIRTH = 0 
const RXN_MIG = 1 
const RXN_DEATH = 2 
const RXN_DYNVAR = 3 

# """
# 	Reaction
#
# Represents a reaction in the structured FGY model.
#
# # Fields
# - `source::String`: Source deme or variable of the reaction
# - `recipient::Union{Nothing,String}`: Recipient deme or variable of the reaction (if applicable)
# - `type::Int`: Type of reaction (RXN_BIRTH, RXN_MIG, RXN_DEATH, or RXN_DYNVAR)
# - `expr::Expr`: Expression defining the reaction rate
# """
mutable struct Reaction
	source::String
	recipient::Union{Nothing,String}
	type::Int
	expr::Expr 
	function Reaction(s::String,r::Union{Nothing,String},t::Int,e::Union{Symbol,Expr})
		@assert t in [ RXN_BIRTH, RXN_MIG, RXN_DEATH, RXN_DYNVAR ]
		o = new()
		o.source = s; o.recipient = r; o.type = t; 
		o.expr = typeof(e)===Symbol ? Expr(e) : e 
		o
	end
end

# """
# 	Reaction(s::String, t::Int, e::Union{Symbol,Expr})
#
# Constructor for Reaction without a recipient (for death and dynamic variable reactions).
#
# # Arguments
# - `s::String`: Source deme or variable
# - `t::Int`: Type of reaction (should be RXN_DEATH or RXN_DYNVAR)
# - `e::Union{Symbol,Expr}`: Expression defining the reaction rate
# """
function Reaction(s::String, t::Int, e::Union{Symbol,Expr})
	@assert t in [ RXN_DEATH, RXN_DYNVAR ]
	Reaction(s, nothing, t
	,  typeof(e)===Symbol ? Expr(e) : e  
	)
end

"""
	ModelFGY

Represents a structured Forward-in-time Genealogy (FGY) model for coalescent simulations.

# Fields
- `modelname::String`: Name of the model
- `birthrxn::Array{Reaction}`: Birth reactions in the model
- `migrationrxn::Array{Reaction}`: Migration reactions in the model
- `deathrxn::Array{Reaction}`: Death reactions in the model
- `nondemerxn::Array{Reaction}`: Non-deme reactions in the model
- `demes::Array{String}`: Names of demes in the model
- `nondemes::Union{Nothing,Array{String}}`: Names of non-deme variables
- `numberdemes::Int`: Number of demes
- `numbernondemes::Int`: Number of non-deme variables
- `initial::Dict{String,Number}`: Initial conditions for variables
- `t0::Float64`: Initial time
- `tfin::Float64`: Final time
- `parameters::Union{Nothing,Dict{String,Number}}`: Model parameters
- `helperexprs::Union{Nothing,Array{Expr}}`: Helper expressions for the model
"""
mutable struct ModelFGY
	modelname::String
	birthrxn::Array{Reaction}
	migrationrxn::Array{Reaction}
	deathrxn::Array{Reaction}
	nondemerxn::Array{Reaction}
	demes::Array{String}
	nondemes::Union{Nothing,Array{String}}
	numberdemes::Int
	numbernondemes::Int 
	initial::Dict{String,Number}
	t0::Float64
	tfin::Float64
	parameters::Union{Nothing,Dict{String,Number}}
	helperexprs::Union{Nothing,Array{Expr}}
end

function Base.show(io::IO, x::ModelFGY)
	print("""
Compartmental model with $(x.numberdemes) demes, 
and $(x.numbernondemes) other dynamic variables.

Dynamic variables: $(x.demes), $(x.nondemes)

Parameters: 
$(DataFrame(Tuple(x.parameters), ["parameter", "value"]))

Initial conditions: $(x.initial)
$(DataFrame(Tuple(x.initial), ["variable", "initial value"]))

Initial time: $(x.t0)

Final time: $(x.tfin)

""")
end
#= 
"""
	ModelFGY(; kwargs...)

Constructor for ModelFGY using keyword arguments.

# Keywords
- `modelname::String`: Name of the model
- `birthrxn::Array{Reaction}`: Birth reactions
- `migrationrxn::Array{Reaction}`: Migration reactions
- `deathrxn::Array{Reaction}`: Death reactions
- `nondemerxn::Array{Reaction}`: Non-deme reactions
- `demes::Array{String}`: Names of demes
- `nondemes::Union{Nothing,Array{String}}`: Names of non-deme variables
- `initial::Dict{String,Float64}`: Initial conditions
- `t0::Float64`: Initial time
- `tfin::Float64`: Final time
- `parameters::Dict{String,Float64}`: Model parameters
- `helperexprs::Union{Nothing,Array{Expr}}`: Helper expressions
"""
function ModelFGY(
		  modelname::String
		, birthrxn::Array{Reaction}
		, migrationrxn::Array{Reaction}
		, deathrxn::Array{Reaction}
		, nondemerxn::Array{Reaction}
		, demes::Array{String}
		, nondemes::Union{Nothing,Array{String}}
		, initial::Dict{String,Float64}
		, t0::Float64
		, tfin::Float64
		, parameters::Dict{String,Float64}
		, helperexprs::Union{Nothing,Array{Expr}}
	)
	ModelFGY(
 		 modelname
		, birthrxn
		, migrationrxn
		, deathrxn
		, nondemerxn
		, demes 
		, nondemes 
		, length(demes) 
		, length(nondemes)
		, initial 
		, t0 
		, tfin 
		, parameters 
		, helperexprs
	)
end =#

"""
	ModelFGY(conffn::String)

Constructor for ModelFGY from a YAML configuration file.

# Arguments
- `conffn::String`: Path to the YAML configuration file

# Returns
- `ModelFGY`: The constructed model
"""
function ModelFGY(conffn::String)
	ModelFGY(confstr = read(conffn, String) )
end 


"""
	ModelFGY(; confstr::String)

Constructor for ModelFGY from a YAML configuration string.

# Arguments
- `confstr::String`: String defining model in YAML format 

# Returns
- `ModelFGY`: The constructed model
"""
function ModelFGY(; confstr::String)
	# conf = YAML.load_file(conffn)
	conf = YAML.load(confstr)
	
	modelname = conf["modelname"] 
	brxns = map(b -> Reaction(b["source"]
		    , b["recipient"]
		    , RXN_BIRTH
		    , Meta.parse(b["rate"]))
	  , conf["births"])
	migrxns = "migrations" ∈ keys(conf) ? map(b -> Reaction(b["source"]
		    , b["recipient"]
		    , RXN_MIG
		    , Meta.parse(b["rate"]))
	  	    , conf["migrations"]) : Array{Reaction}(undef, 0)
	deathrxns = "deaths" ∈ keys(conf) ?  map(b -> Reaction(b["deme"], RXN_DEATH, Meta.parse(b["rate"]))
		, conf["deaths"]) : Array{Reaction}(undef, 0)	
	
	demes = Array{String}(union([rxn.source for rxn in brxns]
		  , [rxn.recipient for rxn in brxns]
	          , [rxn.source for rxn in migrxns]
	          , [rxn.recipient for rxn in migrxns] 
	          , [rxn.source for rxn in deathrxns]
			  ))
	
	numberdemes = length(demes)
	
	dvdict = Dict(zip([x["name"] for x in conf["dynamic_variables"]], conf["dynamic_variables"]))
	dvkeys = keys(dvdict)
	dvvals = [v["initial_value"]  for (k,v) in dvdict]
	initials = Dict(zip(dvkeys, dvvals))
	
	nondemes = [x for x in setdiff(keys(initials), demes)]
	numbernondemes = length(nondemes)
	nondemerxn = Reaction[]
	for k in nondemes
		if "ode" ∉ keys(dvdict[k])
			@warn "Missing ODE definition for dynamic variable $k"
			push!(nondemerxn, Reaction(k, RXN_DYNVAR, :(0+0)))
		else 
			push!(nondemerxn, Reaction(k, RXN_DYNVAR, Meta.parse(dvdict[k]["ode"])))
		end
	end

	t0 = float(conf["time"]["initial"])
	tfin = float(conf["time"]["final"])
	
	parmdict = [(d["name"],d["value"]) for d in conf["parameters"] if d["value"] isa Number] |> Dict{String,Float64}
	strvalparms = [(d["name"], eval(Meta.parse(d["value"]))) for d in conf["parameters"] if d["value"] isa String]
	if length(strvalparms) > 0 
		parmdict = merge(parmdict, strvalparms |> Dict{String,Float64})
	end

	if "helpers" ∈ keys(conf)
		for d in conf["helpers"]
			@assert d["definition"] isa String
		end
	end
	helperexprs = "helpers" ∈ keys(conf) ? [:($(Symbol(d["name"])) = $(Meta.parse(d["definition"]))) for d in conf["helpers"]] :  nothing 
	# ModelFGY( 
	# 	modelname = modelname
	# 	, birthrxn = brxns
	# 	, migrationrxn = migrxns
	# 	, deathrxn = deathrxns
	# 	, nondemerxn = nondemerxn
	# 	, demes = demes 
	# 	, nondemes = nondemes 
	# 	, initial = initials 
	# 	, t0 = t0 
	# 	, tfin = tfin 
	# 	, parameters = parmdict 
	# 	, helperexprs = helperexprs 
	# )
	ModelFGY( 
		modelname
		,  brxns
		,  migrxns
		,  deathrxns
		,  nondemerxn
		,  demes 
		,  nondemes 
		, length(demes) 
		, length(nondemes)
		,  initials 
		,  t0 
		,  tfin 
		,  parmdict 
		,  helperexprs 
	)
end

# """
# 	_derive_deme_ode(deme::String, mod::ModelFGY)
#
# Derive the ODE for a specific deme in the model.
#
# # Arguments
# - `deme::String`: Name of the deme
# - `mod::ModelFGY`: The model
#
# # Returns
# - `Expr`: The derived ODE expression for the deme
# """
function _derive_deme_ode(deme::String, mod::ModelFGY)
	inbrxn = [r.expr for r in mod.birthrxn  if r.recipient==deme]
	inbrxn = length(inbrxn) == 0 ? [:(0)] : inbrxn 
	inmigrxn = [r.expr for r in mod.migrationrxn  if r.recipient==deme]
	inmigrxn = length(inmigrxn) == 0 ? [:(0)] : inmigrxn 
	outmigrxn = [r.expr for r in mod.migrationrxn  if r.source==deme]
	outmigrxn = length(outmigrxn) == 0 ? [:(0)] : outmigrxn 
	deathrxn = [r.expr for r in mod.deathrxn if r.source==deme]
	deathrxn = length(deathrxn) == 0 ? [:(0)] : deathrxn 
	p = Expr(:call,:+, inbrxn..., inmigrxn...) 
	n = Expr(:call,:+, outmigrxn..., deathrxn...)
	Expr(:call, :-, p, n)
end

"""
	solveodes(model::ModelFGY; odemethod = :(Rosenbrock23()), res::Union{Missing,Int64} = missing)

Solve the ordinary differential equations (ODEs) defined by the given model.

# Arguments
- `model::ModelFGY`: The structured FGY model containing the ODEs to solve

# Keywords
- `odemethod = :(Rosenbrock23())`: The ODE solver method to use
- `res::Union{Missing,Int64} = missing`: The number of time points to return in the solution

# Returns
- `ODESolution`: The solution of the ODEs
"""
function solveodes(model::ModelFGY; odemethod = :(Rosenbrock23()), res::Union{Missing,Int64} = missing)
	odees = vcat( 
		map(d ->  _derive_deme_ode(d,model), model.demes) 
		, map(r -> r.expr, model.nondemerxn)
		)
	odeexprs = [:(du[$i] = $ode) for (i,ode) in enumerate(odees)]
	odeexpr = Expr(:block, odeexprs...)

	assexprs = [:($(Symbol(deme)) = u[$(i)]) for (i,deme) in enumerate(vcat(model.demes,model.nondemes))]
	assexpr = Expr(:block, assexprs...)

	helperexpr = isnothing(model.helperexprs) ? :() : Expr(:block, model.helperexprs...)
	
	paex =  [(:($(Symbol(k)) = $v)) for (k,v) in model.parameters]
	for ex in paex 
		eval(ex) 
	end

	eval(quote 
		function mododes!(du, u, p, t)
			$assexpr
			$helperexpr
			$odeexpr
		end
	end)
	
	initial_cond = [model.initial[k] for k in vcat(model.demes, model.nondemes)]
	trange = (model.t0, model.tfin)
	
	prmod = ODEProblem(mododes! 
		, initial_cond 
		, trange
	)
	integ = eval(odemethod)
	tstops = ismissing(res) ? [] : collect(range(model.t0, model.tfin, length=res))

	smod = solve(prmod, integ; tstops = tstops)
	smod 
end
