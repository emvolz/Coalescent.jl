#= 
Structured FGY model 
=#

# using .Coalescent
using Debugger 
using YAML 
using DataFrames
using OrdinaryDiffEq
using Plots

const RXN_BIRTH = 0 
const RXN_MIG = 1 
const RXN_DEATH = 2 
const RXN_DYNVAR = 3 

mutable struct Reaction
	source::String
	recipient::Union{Nothing,String}
	type::Int
	expr::Expr 
	function Reaction(s::String,r::Union{Nothing,String},t::Int,e::Expr)
		@assert t in [ RXN_BIRTH, RXN_MIG, RXN_DEATH, RXN_DYNVAR ]
		o = new()
		o.source = s; o.recipient = r; o.type = t; o.expr = e
		o
	end
end
function Reaction( s::String, t::Int, e::Expr )
	@assert t in [ RXN_DEATH, RXN_DYNVAR ]
	Reaction( s, nothing, t,  e )
end

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
	parameters::Dict{String,Number}
	helperexprs::Array{Expr}
end
function Base.show( io::IO, x::ModelFGY)
	print("""
Compartmental model with $(x.numberdemes) demes, 
and $(x.numbernondemes) other dynamic variables.

Dynamic variables: $(x.demes), $(x.nondemes)

Parameters: 
$(DataFrame( Tuple( x.parameters ), ["parameter", "value"] ))

Initial conditions: $(x.initial)
$(DataFrame( Tuple( x.initial ), ["variable", "initial value"] ))

Initial time: $(x.t0)

Final time: $(x.tfin)

""")
end

function ModelFGY(;
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
		, helperexprs::Array{Expr}
	)
	ModelFGY(
 		 modelname
		, birthrxn
		, migrationrxn
		, deathrxn
		, nondemerxn
		, demes 
		, nondemes 
		, length( demes ) 
		, length( nondemes )
		, initial 
		, t0 
		, tfin 
		, parameters 
		, helperexprs
	)
end
function ModelFGY(conffn::String)
	conf = YAML.load_file( conffn )
	
	modelname = conf["modelname"] 
	brxns = map( b -> Reaction(b["source"]
		    , b["recipient"]
		    , RXN_BIRTH
		    , Meta.parse(b["rate"]))
	  , conf["births"] )
	migrxns = map( b -> Reaction(b["source"]
		    , b["recipient"]
		    , RXN_MIG
		    , Meta.parse(b["rate"]))
	  , conf["migrations"] )
	# TODO  migrations should be optional
	# deaths should be optional ( but maybe raise warnings if missing) 
	deathrxns = map( b -> Reaction(b["deme"], RXN_DEATH, Meta.parse(b["rate"]) )
	  , conf["deaths"] )
	
	demes = union( [ rxn.source for rxn in brxns ]
		  , [ rxn.recipient for rxn in brxns ]
	          , [ rxn.source for rxn in migrxns ]
	          , [ rxn.recipient for rxn in migrxns ] 
	          , [ rxn.source for rxn in deathrxns ]
	)
	#|> setdiff( nothing )
	
	numberdemes = length( demes )
	
	dvdict = Dict( zip( [x["name"] for x in conf["dynamic_variables"]], conf["dynamic_variables"] ) )
	dvkeys = keys(dvdict)
	dvvals = [ v["initial_value"]  for (k,v) in dvdict ]
	initials = Dict( zip( dvkeys, dvvals ))
	
	nondemes = [ x for x in setdiff( keys( initials ), demes ) ]
	numbernondemes = length( nondemes )
	nondemerxn = [ Reaction(k, RXN_DYNVAR, Meta.parse(dvdict[k]["ode"])) for k in nondemes]  

	t0 = float( conf[ "time"]["initial"] )
	tfin = float( conf[ "time"]["final"] )
	
	parmdict = [ (d["name"],d["value"]) for d in conf["parameters"] if d["value"] isa Number ] |> Dict
	parmdict = merge( parmdict, [ (d["name"], eval(Meta.parse(d["value"]))) for d in conf["parameters"] if d["value"] isa String ] |> Dict )

	if "helpers" ∈ keys(conf)
		for d in conf["helpers"]
			@assert d["definition"] isa String
		end
	end
	helperexprs = "helpers" ∈ keys(conf) ? [ :($(Symbol(d["name"])) = $(Meta.parse(d["definition"])) ) for d in conf[ "helpers" ] ] : []  

	ModelFGY( 
		modelname = modelname
		, birthrxn = brxns
		, migrationrxn = migrxns
		, deathrxn = deathrxns
		, nondemerxn = nondemerxn
		, demes = demes 
		, nondemes = nondemes 
		, initial = initials 
		, t0 = t0 
		, tfin = tfin 
		, parameters = parmdict 
		, helperexprs = helperexprs 
	)
end

function _derive_deme_ode( deme::String, mod::ModelFGY )
	inbrxn = [ r.expr for r in mod.birthrxn  if r.recipient==deme ]
	inmigrxn = [ r.expr for r in mod.migrationrxn  if r.recipient==deme ]
	outmigrxn = [ r.expr for r in mod.migrationrxn  if r.source==deme ]
	deathrxn = [ r.expr for r in mod.deathrxn if r.source==deme ]
	p = Expr(:call,:+, inbrxn..., inmigrxn...) 
	n = Expr(:call,:+, outmigrxn..., deathrxn...)
	Expr(:call, :-, p, n )
end

function solveodes(model::ModelFGY)
	odees = vcat( 
		map( d ->  _derive_deme_ode(d,model), model.demes ) 
		, map( r -> r.expr, model.nondemerxn )
		)
	odeexprs = [ :(du[$i] = $ode) for (i,ode) in enumerate(odees)]
	odeexpr = Expr(:block, odeexprs... )

	assexprs = [ :( $(Symbol(deme)) = u[$(i)]) for (i,deme) in enumerate( vcat(model.demes,model.nondemes) ) ]
	assexpr = Expr( :block, assexprs... )

	helperexpr = Expr(:block, model.helperexprs... )
	
	# TODO move inside mododes to prevent name conflict 
	paex =  [ ( :($(Symbol(k)) = $v) ) for (k,v) in model.parameters ]
	for ex in paex 
		# @eval ex 
		eval(ex) 
	end

	eval( quote 
		function mododes!( du, u, p, t )
			$assexpr
			$helperexpr
			$odeexpr
		end
	end )
	
	initial_cond = [ model.initial[k] for k in vcat(model.demes, model.nondemes) ]
	trange = ( model.t0, model.tfin )
	
	prmod = ODEProblem( mododes! 
		, initial_cond 
		, trange
	)
	integ = Rosenbrock23() 
	smod = solve( prmod, integ )
	smod 
end

