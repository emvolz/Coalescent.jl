using YAML 
using DataFrames
using OrdinaryDiffEq
using Statistics
using Random
using Distributions
using Interpolations
using Plots
using StatsBase
import MacroTools
using JumpProcesses

"""
    SimTree(events::Array{Event}, model::ModelFGY; computedescendants = false)::SimTree

Construct a SimTree from an array of Events and a ModelFGY.

# Arguments
- `events::Array{Event}`: Array of coalescent events
- `model::ModelFGY`: The structured FGY model

# Keywords
- `computedescendants::Bool = false`: Whether to compute descendants for each node

# Returns
- `SimTree`: The constructed simulated tree
"""
function SimTree(events::Array{Event}, model::ModelFGY; computedescendants = false)::SimTree
	
	ix = sortperm( [ e.height for e in events ] )
	events = events[ ix ]
	n = sum( [ e.type == SAMPLE for e in events ] )
	nNode = n - 1 
	nCo = sum( [e.type == COALESCENT for e in events ] )
	mrcaheight = events[end].height 
	epsilonheight = mrcaheight - events[end-1].height 
	cotoadd = (nNode - nCo )
	if cotoadd > 0 
		@warn "Coalescent process did not reach a common ancestor. Adding $(cotoadd) nodes. "
	end
	for i in 1:cotoadd # add coalescent events if tree is not complete (single mrca)
		push!( events, Event(COALESCENT, mrcaheight + epsilonheight*i) )
	end
	
	deme2n = Dict(zip(model.demes,
		   [ sum(map(e -> e.type==SAMPLE && e.source==deme, events )) 
			for deme in model.demes ]
	))
	
	tiplabs = [ "t.$(e.source).$(i)" for (i,e) in enumerate(events) if e.type == SAMPLE ]
	shs =  [ e.height for e in events if e.type == SAMPLE ]
	nodes = range( 1, n + nNode ) |> collect 
	

	nedges = n + nNode - 1
	parent = Array{Union{Nothing,Int}}(nothing, nedges )
	child = Array{Union{Nothing,Int}}(nothing, nedges )
	daughters = Vector{Tuple{Int, Int, Int, Int }}(undef, nNode  )
	edgelength = Array{Union{Nothing,Float64}}(nothing, nedges )
	heights = Array{Union{Nothing,Float64}}( nothing, n + nNode )
	demes = Array{String}(undef, n + nNode )
	
	A = 0
	ie = 1
	
	extant = Dict( zip( model.demes, [ Array{Int}(undef, 0) for _ in model.demes ] ))
	for v in values(extant)
		sizehint!( v, n + nNode )
	end

	a = n + nNode # internal node counter, goes down 
	su = 1 # tip counter, goes up 

	desc = computedescendants ?  falses( nNode, n ) : nothing 
	
	iu = 0 
	iv = 0
	for e in events
		if e.type == SAMPLE
			A+=1 
			heights[su] = e.height 
			push!( extant[e.source], su )
			demes[su] = e.source 
			su += 1
		elseif e.type == COALESCENT
			heights[a] = e.height

			if isnothing( e.source) # random coalesce at mrcaheight to get binary tree 
				# Adict = Dict( zip( model.demes,  map(x -> length(extant[x]), model.demes) ))
				Avec =  map(x -> length(extant[x]), model.demes)
				demeu = sample( model.demes, Weights(Avec), 1 )[1]
				e.source = demeu 
				Avec[ model.demes .== demeu] .-= 1 
				# Adict[ demeu ] -= 1 
				demev = sample( model.demes, Weights(Avec), 1 )[1]
				e.sink = demev 
				Avec[ model.demes .== demev] .-= 1 
			end

			e0 = extant[ e.source ]
			e1 = extant[ e.sink ]

			iu = sample(1:length(e0))
			if e.source == e.sink
				iv = sample( setdiff(1:length(e1), iu) )
			else
				iv = sample(1:length(e1))	
			end
			u = extant[e.source][iu]
			v = extant[e.sink][iv] 
			if e.source == e.sink 
				deleteat!( extant[e.source], sort!([iu,iv]) )
			else 
				deleteat!( extant[e.source], iu )
				deleteat!( extant[e.sink], iv )
			end

			push!(extant[e.source], a )
			demes[a] = e.source 

			A-=1
			parent[ie] = a 
			child[ie] = u 
			edgelength[ie] = heights[a] - heights[u] 
			ie +=1 
			parent[ie] = a 
			child[ie] = v 
			edgelength[ie] = heights[a] - heights[v]
			ie += 1 
			daughters[ a - n ] = (u,v, ie-2, ie-1) # edgelength[ie-2], edgelength[ie-1]) 
			if computedescendants
				if u <= n 
					desc[a-n,u] = true 
				else 
					desc[a-n,:] .|= desc[u-n,:] 
				end
				if v <= n 
					desc[a-n,v] = true 
				else 
					desc[a-n,:] .|= desc[v-n,:] 
				end
			end
			a -= 1
		elseif e.type == MIGRATION
			e1 = extant[ e.sink ]
			iv = sample(1:length(e1))
			v = e1[iv] 
			deleteat!( e1, iv )
			push!( extant[e.source], v )
		else
			Error( "Not Implemented" )
		end
	end
	
	SimTree( parent, child, n, nNode, edgelength, heights, tiplabs, shs, desc, daughters, demes )
end

"""
    SimTree(model::ModelFGY, sample::SampleConfiguration; computedescendants = false)

Simulate a coalescent tree based on a given model and sampling configuration.

# Arguments
- `model::ModelFGY`: The structured FGY model to simulate
- `sample::SampleConfiguration`: The sampling configuration for the simulation

# Keywords
- `computedescendants::Bool = false`: Whether to compute descendants for each node

# Returns
- `SimTree`: A simulated coalescent tree
"""
function SimTree( model::ModelFGY, sample::SampleConfiguration; computedescendants = false  )
	_sim_markov( model
	  , [ x[2] for x in sample.sconf ]
	  , [ x[1] for x in sample.sconf ]
	  , computedescendants 
	)
end



# """
#     _sim_markov(model::ModelFGY, sampletimes::Array{Float64}, samplestates::Array{String}, computedescendants::Bool; odemethod = :(AutoTsit5(Rosenbrock23())), ytol = 1e-6)
#
# Simulate a coalescent tree using a Markovian coalescent algorithm for a structured population model.
#
# # Arguments
# - `model::ModelFGY`: The structured FGY model to simulate
# - `sampletimes::Array{Float64}`: Array of sample times
# - `samplestates::Array{String}`: Array of sample states (demes)
# - `computedescendants::Bool`: Whether to compute descendants for each node
#
# # Keywords
# - `odemethod = :(AutoTsit5(Rosenbrock23()))`: The ODE solver method to use
# - `ytol::Float64 = 1e-6`: Tolerance for numerical stability in rate calculations
#
# # Returns
# - `SimTree`: The simulated coalescent tree
# """
function _sim_markov( model::ModelFGY
		     , sampletimes::Array{Float64}
		     , samplestates::Array{String}
		     #, tmrcaguess::Float64
		     #, p... 
		     , computedescendants::Bool 
		     ; odemethod = :(AutoTsit5(Rosenbrock23()))
		     , ytol = 1e-6
	)
	n = length( sampletimes )
	@assert n > 1
	# @assert tmrcaguess > 0. 
	# @assert diff( collect( extrema( sampletimes)))[1] < tmrcaguess 
	@assert length( sampletimes ) == length( samplestates )
	
	# sort!(sampletimes)
	mst = maximum( sampletimes )
	ix = sortperm( sampletimes, rev = true )
	sampleheights = mst .- sampletimes[ix]  
	samplestates = samplestates[ix] 
	# samplesadded = fill( false, n )
	
	ushs = unique( sampleheights ) |> sort!

	deme2shsdict = map(deme -> countmap( sampleheights[ samplestates.==deme ] ), model.demes )
	currentsampleheight = 0.0
	ixsampleheight = 1

	events = Array{Event}(undef, 0)
	sizehint!(events, n*10)

	# solve model 
	msol = solveodes( model )

	((msol.t[end] >= mst) || (msol.t[end] ≈ mst)) ||  throw(error( "Some sample times exceed final simulation time point." ))


	
	# interpolators 
	_fnmsolu(i) = [ x[i] for x in msol.u ]
	interpdict = Dict( zip( 
		vcat( model.demes , model.nondemes ) # columns of msol.u are in this order
		, [ linear_interpolation( reverse(mst.-msol.t), reverse(_fnmsolu(i)) )
			for (i,k) in enumerate( vcat( model.demes , model.nondemes ) ) ]
	))

	# sampadds = 0 
	cou = rand() 
	sampcount = 0 
	cocount = 0 
	coheights = Array{Float64}(undef,0); sizehint!(coheights, n-1)
	ndemes = length( model.demes )
	A = Dict( zip( model.demes, fill(0., ndemes ) ))

	assexprs = [ :( $(Symbol(v)) = interpdict[$v](t)) for (i,v) in enumerate( vcat(model.demes, model.nondemes) ) ]
	assexpr = Expr( :block, assexprs... )
	
	mparameters = isnothing( model.parameters ) ? Dict{String,Float64}() : model.parameters
	mparameters["mst"] = mst # NOTE copying here so the variable is accessible inside ODE expressions 
	paex =  Expr( :block, [ ( :($(Symbol(k)) = $v) ) for (k,v) in mparameters]... )
	helperexpr = isnothing(model.helperexprs) ? :() :  Expr( :block, model.helperexprs... )
	Aex = Expr( :block, [ ( :($(Symbol("A_"*k)) = A[$k]) ) for k in keys(A)]... )

	"Replace occurrences of the time variable 't' with '(mst-t)' in the given expression."
	function _sanitize_time_variable( expr )
		MacroTools.postwalk( expr ) do x 
			x == :t ? :(mst-t) : x 
		end
	end

	"Generate the expression for coalescence rate."
	function _corateex( rxn ; ytol = ytol )
		@assert rxn.type == RXN_BIRTH
		@assert !isnothing( rxn.recipient )
		src = rxn.source 
		rec = rxn.recipient

		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)
		Yₛ = :( max( $Aₛ*$ytol,   $(Symbol(src))))
		Yᵣ = :( max( $Aᵣ*$ytol,   $(Symbol(rec))))

		if src == rec 
			aex = :( $Aₛ * max(0.0,$Aᵣ-1.0) / ($Yₛ * $Yᵣ) )
		else
			aex = :( $Aₛ * $Aᵣ / ($Yₛ * $Yᵣ) )
		end
		Expr( :call, :*, aex,  _sanitize_time_variable( rxn.expr ) )
	end

	"Generate the expression for birth-migration rate."
	function _birthmigrateex( rxn ; ytol = ytol  )
		# forward time transmission s -> r, reverse time migration r -> s
		@assert rxn.type == RXN_BIRTH
		@assert !isnothing( rxn.recipient )
		src = rxn.source 
		rec = rxn.recipient

		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)
		Yₛ = :( max( $Aₛ*$ytol,   $(Symbol(src))))
		Yᵣ = :( max( $Aᵣ*$ytol,   $(Symbol(rec))))

		pex = :( clamp( ($Yₛ - $Aₛ) / $Yₛ , 0., 1. ) )
		aex = Expr( :call, :*, :( $Aᵣ / $Yᵣ ), pex )

		Expr( :call, :*, aex,  _sanitize_time_variable( rxn.expr )  )
	end

	"Generate the expression for migration rate."
	function _migrateex( rxn ; ytol = ytol )
		# reverse time migration r -> s
		@assert rxn.type == RXN_MIG
		@assert rxn.source != rxn.recipient 
		src = rxn.source 
		rec = rxn.recipient
		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)
		# Yₛ = :( max( 0.0, $(Symbol(src))))
		# Yᵣ = :( max(0.0, $(Symbol(rec))))
		Yₛ = :( max( $Aₛ*$ytol,  $(Symbol(src))))
		Yᵣ = :( max( $Aᵣ*$ytol,   $(Symbol(rec))))

		aex = :( $Aᵣ / $Yᵣ ) 
		
		Expr( :call, :*, aex,  _sanitize_time_variable( rxn.expr ) )
	end

	birthmigrxns = [ r for r in model.birthrxn if r.source != r.recipient ]
	coexprs = map( r -> _corateex(r), model.birthrxn )
	birthmigexprs = map( r -> _birthmigrateex(r), birthmigrxns )
	migexprs = map( r -> _migrateex(r), model.migrationrxn )
	eventexprs = vcat( coexprs , birthmigexprs , migexprs )
	eventtypes = vcat( fill(COALESCENT,length(coexprs)) 
		, fill( MIGRATION, length(birthmigexprs)) 
		, fill( MIGRATION, length(migexprs))
	)
	eventrxns = vcat( model.birthrxn ,  birthmigrxns,  model.migrationrxn )
	nevents = length(eventexprs)
	eventrates = fill( 0.0, nevents)
	eventexprs1 = [ :(eventrates[$i] = $ex) for (i,ex) in enumerate(eventexprs) ]
	eventexprblock = Expr(:block, eventexprs1... )

	eval( quote
	function fneventrates!(eventrates, t, A, interpdict) 
		$Aex 
		$assexpr 
		$paex
		$helperexpr
		$eventexprblock
		eventrates[ isnan.(eventrates) ] .= 0.0  
		nothing
	end
	end )
	
	function coodes!(du, u, p, t)
		# fneventrates!(eventrates, t, A, interpdict) # A, interpdict 
		Base.invokelatest( fneventrates!, eventrates, t, A, interpdict )
		du[1] = max(0., sum(eventrates) ); 
	end

	"Trigger callback." 
	function sampcondition(u,t,integrator )::Bool
		t == ushs[ ixsampleheight ] 
	end
	"Update state based on sampling." 
	function sampaffect!(integrator)
		currentsampleheight = ushs[ ixsampleheight ] 
		for (ideme,deme) in enumerate( model.demes )
			shsdict = deme2shsdict[ideme] 
			if currentsampleheight in keys(shsdict)
				newsamps  = shsdict[currentsampleheight]
				A[deme] += newsamps
				append!(events, 
					fill( Event(SAMPLE
			  	  	  , currentsampleheight
			  	  	  , deme 
			  	  	  , deme 
					), newsamps )
				)
			end
		end
		ixsampleheight = min( length( ushs ), ixsampleheight + 1 ) 
	end
	sampcb = DiscreteCallback( sampcondition, sampaffect! )
	
	"Triggers callback. Check if an event should occur."
	function eventcondition(u,t,integrator)::Real
		exp(-u[1]) - cou # note: not integrator.u[1]
	end

	"Handle the occurrence of an event."
	function eventaffect!(integrator)
		# fneventrates!(eventrates, integrator.t, A, interpdict) # A, interpdict 
		Base.invokelatest( fneventrates!, eventrates, integrator.t, A, interpdict )
		we = sample( Weights(eventrates) )
		rxntype = eventtypes[ we ]
		rxn = eventrxns[we]
		if rxntype == COALESCENT
			A[ rxn.recipient ] -= 1
			push!(events, Event(COALESCENT
			  , integrator.t 
			  , rxn.source
			  , rxn.recipient )
			)
		elseif rxntype == MIGRATION
			A[ rxn.recipient ] -= 1
			A[ rxn.source ] += 1
			push!(events, Event(MIGRATION
			  , integrator.t 
			  , rxn.source
			  , rxn.recipient )
			)
		end
		integrator.u[1] = 0. 
		cou = rand()
	end
	eventcb = ContinuousCallback(eventcondition, eventaffect! )

	cbs = CallbackSet( sampcb, eventcb )
	
	# integ = Rosenbrock23()
	# integ = AutoTsit5(Rosenbrock23())
	# integ = Tsit5()
	integ = eval( odemethod )
	sampaffect!( integ ) #initial sample 
	pr = ODEProblem( coodes!
		 , [0.0]
		 , ( 0., mst - model.t0)
		)
	s = solve( pr, integ 
	, callback = cbs, tstops = ushs[ushs.>0.0]  )
	
	SimTree( events, model; computedescendants = computedescendants)
end

