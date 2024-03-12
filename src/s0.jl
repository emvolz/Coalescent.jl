using Coalescent
using YAML 
using DataFrames
using OrdinaryDiffEq
using Statistics
using Random
using Distributions
using RCall
using Interpolations
using Plots
using StatsBase
using Debugger

function SimTree(events::Array{Event}, model::ModelFGY)::SimTree
	ix = sortperm( [ e.height for e in events ] )
	events = events[ ix ]
	n = sum( [ e.type == SAMPLE for e in events ] )
	# coes = [ e for e in events if e.type == COALESCENT ]
	nNode = sum( [e.type == COALESCENT for e in events ] )
	deme2n = Dict(zip(model.demes,
		   [ sum(map(e -> e.type==SAMPLE && e.source==deme, events )) 
			for deme in model.demes ]
	))
	
	tiplabs = [ "t.$(e.source).$(i)" for (i,e) in enumerate(events) if e.type == SAMPLE ]
	nodes = range( 1, n + nNode ) |> collect 
	
	nedges = n + nNode - 1
	parent = Array{Union{Nothing,Int}}(nothing, nedges )
	child = Array{Union{Nothing,Int}}(nothing, nedges )
	edgelength = Array{Union{Nothing,Float64}}(nothing, nedges )
	heights = Array{Union{Nothing,Float64}}( nothing, n + nNode )
	
	A = 0
	ie = 1
	
	extant = Dict( zip( model.demes, [ Array{Int}(undef, 0) for _ in model.demes ] ))
	for v in values(extant)
		sizehint!( v, n + nNode )
	end

	a = n + nNode # internal node counter, goes down 
	su = 1 # tip counter, goes up 
	
	iu = 0 
	iv = 0
	for e in events
		if e.type == SAMPLE
			A+=1 
			heights[su] = e.height 
			push!( extant[e.source], su )
			su += 1
		elseif e.type == COALESCENT
			heights[a] = e.height
			e0 = extant[ e.source ]
			e1 = extant[ e.sink ]
if length(e1)==0
				@bp
end

			iu = sample(1:length(e0))
			if e.source == e.sink
				iv = sample( setdiff(1:length(e1), iu) )
			else
				iv = sample(1:length(e1))	
			end
			u = extant[e.source][iu]
			v = extant[e.sink][iv] 
			deleteat!( extant[e.source], iu )
			deleteat!( extant[e.sink], iv )
			push!(extant[e.source], a )

			A-=1
			parent[ie] = a 
			child[ie] = u 
			edgelength[ie] = heights[a] - heights[u] 
			ie +=1 
			parent[ie] = a 
			child[ie] = v 
			edgelength[ie] = heights[a] - heights[v]
			ie += 1 
			a -= 1
		elseif e.type == MIGRATION
			e1 = extant[ e.sink ]
if length(e1)==0
				@bp
end
			iv = sample(1:length(e1))
			v = e1[iv] 
			deleteat!( e1, iv )
			push!( extant[e.source], v )
		else
			Error( "Not Implemented" )
		end
	end
	
	SimTree( parent, child, n, nNode, edgelength, heights, tiplabs )
end

function _sim_markov( model::ModelFGY
		     , sampletimes::Array{Float64}
		     , samplestates::Array{String}
		     #, tmrcaguess::Float64
		     #, p... 
	)
	print( """Simulating coalescent, sample size = $(length(sampletimes))
Markovian coalescent algorithm
User-specified model
""" * string(model) ) #$(string(model))
	n = length( sampletimes )
	@assert n > 1
	# @assert tmrcaguess > 0. 
	# @assert diff( collect( extrema( sampletimes)))[1] < tmrcaguess 
	@assert length( sampletimes ) == length( samplestates )
	sort!(sampletimes)
	mst = maximum( sampletimes )
	sampleheights = reverse(  mst .- sampletimes ) 
	# samplesadded = fill( false, n )
	
	ushs = unique( sampleheights ) |> sort!

	deme2shsdict = map(deme -> countmap( sampleheights[ samplestates.==deme ] ), model.demes )
	currentsampleheight = 0.0
	ixsampleheight = 1

	events = Array{Event}(undef, 0)
	sizehint!(events, n*10)

	# solve model 
	msol = solveodes( model )
	@assert msol.t[end] >= maximum(sampletimes)
	
	# interpolators 
	_fnmsolu(i) = [ x[i] for x in msol.u ]
	interpdict = Dict( zip( 
		model.demes ∪ model.nondemes # columns of msol.u are in this order
		, [ linear_interpolation( reverse(mst.-msol.t), reverse(_fnmsolu(i)) )
			for (i,k) in enumerate( model.demes ∪ model.nondemes ) ]
	))

	# sampadds = 0 
	cou = rand() 
	sampcount = 0 
	cocount = 0 
	coheights = Array{Float64}(undef,0); sizehint!(coheights, n-1)
	ndemes = length( model.demes )
	A = Dict( zip( model.demes, fill(0., ndemes ) ))

	# initialise A and add initial samples 
	for (ideme,deme) in enumerate( model.demes )
		shsdict = deme2shsdict[ideme] 
		if currentsampleheight in keys(shsdict)
			newsamps = shsdict[currentsampleheight]
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
	ixsampleheight += 1 
	
	assexprs = [ :( $(Symbol(v)) = interpdict[$v](t)) for (i,v) in enumerate( model.demes ∪ model.nondemes ) ]
	assexpr = Expr( :block, assexprs... )
	
	paex =  Expr( :block, [ ( :($(Symbol(k)) = $v) ) for (k,v) in model.parameters ]... )
	helperexpr = Expr( :block, model.helperexprs... )
	Aex = Expr( :block, [ ( :($(Symbol("A_"*k)) = A[$k]) ) for k in keys(A)]... )

	function _corateex( rxn )
		@assert rxn.type == RXN_BIRTH
		@assert !isnothing( rxn.recipient )
		src = rxn.source 
		rec = rxn.recipient

		Yₛ = :( max(0.0, $(Symbol(src))))
		Yᵣ = :( max(0.0, $(Symbol(rec))))
		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)

		if src == rec 
			aex = :( $Aₛ * ($Aᵣ-1.0) / ($Yₛ * $Yᵣ) )
		else
			aex = :( $Aₛ * $Aᵣ / ($Yₛ * $Yᵣ) )
		end
		Expr( :call, :*, aex,  rxn.expr )
	end
	function _birthmigrateex( rxn )
		# forward time transmission s -> r, reverse time migration r -> s
		@assert rxn.type == RXN_BIRTH
		@assert !isnothing( rxn.recipient )
		@assert rxn.source != rxn.recipient
		src = rxn.source 
		rec = rxn.recipient

		Yₛ = :( max( 0.0,  $(Symbol(src))))
		Yᵣ = :( max(0.0, $(Symbol(rec))))
		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)

		pex = :( clamp( ($Yₛ - $Aₛ) / $Yₛ , 0., 1. ) )
		aex = Expr( :call, :*, :( $Aᵣ / $Yᵣ ), pex )

		Expr( :call, :*, aex,  rxn.expr )
	end
	function _migrateex( rxn )
		# reverse time migration r -> s
		@assert rxn.type == RXN_MIG
		@assert rxn.source != rxn.recipient 
		src = rxn.source 
		rec = rxn.recipient
		Yₛ = :( max( 0.0, $(Symbol(src))))
		Yᵣ = :( max(0.0, $(Symbol(rec))))
		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)
		aex = :( $Aᵣ / $Yᵣ ) 
		
		Expr( :call, :*, aex,  rxn.expr )
	end

	coexprs = map( r -> _corateex(r), model.birthrxn )
	birthmigexprs = map( r -> _birthmigrateex(r), model.birthrxn )
	migexprs = map( r -> _migrateex(r), model.migrationrxn )
	eventexprs = vcat( coexprs , birthmigexprs , migexprs )
	eventtypes = vcat( fill(COALESCENT,length(coexprs)) 
		, fill( MIGRATION, length(birthmigexprs)) 
		, fill( MIGRATION, length(migexprs))
	)
	eventrxns = vcat( model.birthrxn ,  model.birthrxn ,  model.migrationrxn )
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
		nothing
	end
	end )

	function coodes!(du, u, p, t)
		fneventrates!(eventrates, t, A, interpdict) # A, interpdict 
		du[1] = max(0., sum(eventrates) ); 
	end

	function sampcondition(u,t,integrator )::Bool
		# t>0 && t in sampleheights
		t == ushs[ ixsampleheight ] 
	end
	function sampaffect!(integrator)
		# currentsampleheight = integrator.t # NOTE integrator should be stopped on each sample height by solver (tstops). 
		# for deme in model.demes 
		# 	newsamps = sum( (samplestates .== deme) .& (sampleheights .== currentsampleheight) )
		# 	A[deme] += newsamps 
		# 	append!(events, 
		# 		fill( Event(SAMPLE
		# 	  	  , integrator.t 
		# 	  	  , deme 
		# 	  	  , deme 
		# 		), newsamps )
		# 	)
		# end
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
		ixsampleheight = min( length( sampleheights ), ixsampleheight + 1 ) 
	end
	sampcb = DiscreteCallback( sampcondition, sampaffect! )
	
	function eventcondition(u,t,integrator)::Real
		#u[1] - cou
		exp(-u[1]) - cou # note: not integrator.u[1]
	end
	function eventaffect!(integrator)
		fneventrates!(eventrates, integrator.t, A, interpdict) # A, interpdict 
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
# println("############")
# 		println(integrator.t)
# println( DataFrame( hcat( eventtypes, eventrates ), :auto ))
# println( DataFrame( A ) )
		integrator.u[1] = 0. 
		cou = rand()
	end
	eventcb = ContinuousCallback(eventcondition, eventaffect! )

	cbs = CallbackSet( sampcb, eventcb )
	
	# integ = Rosenbrock23()
	integ = AutoTsit5(Rosenbrock23())
	# integ = Tsit5()
	pr = ODEProblem( coodes!
		 , [0.0]
		 #, [ 1.0 ]
		 , ( 0., mst - model.t0)
		)
	s = solve( pr, integ 
	   , callback = cbs, tstops = ushs[ ushs .> 0. ] )
	
# for e in events 
# 		println( e  )
#end
#
	SimTree( events, model )
end

