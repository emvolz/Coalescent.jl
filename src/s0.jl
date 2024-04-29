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
import MacroTools
using JumpProcesses

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
			if e.source == e.sink 
				deleteat!( extant[e.source], sort!([iu,iv]) )
			else 
				deleteat!( extant[e.source], iu )
				deleteat!( extant[e.sink], iv )
			end
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

function SimTree( model::ModelFGY, sample::SampleConfiguration )
println("££££££££££££ sim_markov2 ££££££££££££££")
	_sim_markov( model
	  , [ x[2] for x in sample.sconf ]
	  , [ x[1] for x in sample.sconf ]
	)
end


function _sim_markov2( model::ModelFGY
		     , sampletimes::Array{Float64}
		     , samplestates::Array{String}
		     #, tmrcaguess::Float64
		     #, p... 
	)
	println( """Simulating coalescent, sample size = $(length(sampletimes))
££££££££££££ sim_markov2 ££££££££££££££
Markovian coalescent algorithm
User-specified model
""" * string(model) ) #$(string(model))
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
	mindeltat = minimum(diff(msol.t))
	mediandeltat = median(diff(msol.t))

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
	All = countmap( samplestates )
	All2 = Dict( [(k,n) for k in model.demes ] ) # n individuals in every deme, for calculating rate bounds 

	assexprs = [ :( $(Symbol(v)) = interpdict[$v](t)) for (i,v) in enumerate( vcat(model.demes, model.nondemes) ) ]
	assexpr = Expr( :block, assexprs... )
	
	mparameters = model.parameters
	mparameters["mst"] = mst # NOTE copying here so the variable is accessible inside ODE expressions 
	paex =  Expr( :block, [ ( :($(Symbol(k)) = $v) ) for (k,v) in mparameters]... )
	helperexpr = Expr( :block, model.helperexprs... )
	Aex = Expr( :block, [ ( :($(Symbol("A_"*k)) = A[$k]) ) for k in keys(A)]... )

	function _sanitize_time_variable( expr )
		#= 
		if time (t) appears in expr, replace with (mst-t) for coalescent equations
				example input:
				:(if t > 80.6
					mu * hostB3
				else
					0.0
				end)
		=# 
		MacroTools.postwalk( expr ) do x 
			x == :t ? :(mst-t) : x 
		end
	end
	function _corateex( rxn )
		@assert rxn.type == RXN_BIRTH
		@assert !isnothing( rxn.recipient )
		src = rxn.source 
		rec = rxn.recipient

		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)
		Yₛ = :( max( $Aₛ,   $(Symbol(src))))
		Yᵣ = :( max( $Aᵣ,   $(Symbol(rec))))

		if src == rec 
			aex = :( $Aₛ * max(0.0,$Aᵣ-1.0) / ($Yₛ * $Yᵣ) )
		else
			aex = :( $Aₛ * $Aᵣ / ($Yₛ * $Yᵣ) )
		end
		Expr( :call, :*, aex,  _sanitize_time_variable( rxn.expr ) )
	end
	function _birthmigrateex( rxn )
		# forward time transmission s -> r, reverse time migration r -> s
		@assert rxn.type == RXN_BIRTH
		@assert !isnothing( rxn.recipient )
		src = rxn.source 
		rec = rxn.recipient

		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)
		Yₛ = :( max( $Aₛ, $(Symbol(src))))
		Yᵣ = :( max( $Aᵣ, $(Symbol(rec))))

		pex = :( clamp( ($Yₛ - $Aₛ) / $Yₛ , 0., 1. ) )
		aex = Expr( :call, :*, :( $Aᵣ / $Yᵣ ), pex )

		Expr( :call, :*, aex,  _sanitize_time_variable( rxn.expr )  )
	end
	function _migrateex( rxn )
		# reverse time migration r -> s
		@assert rxn.type == RXN_MIG
		@assert rxn.source != rxn.recipient 
		src = rxn.source 
		rec = rxn.recipient
		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)
		# Yₛ = :( max( 0.0, $(Symbol(src))))
		# Yᵣ = :( max( 0.0, $(Symbol(rec))))
		Yₛ = :( max( $Aₛ, $(Symbol(src))))
		Yᵣ = :( max( $Aᵣ, $(Symbol(rec))))
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
	function fneventrates( t, A, interpdict, nevents) 
		$Aex 
		$assexpr 
		$paex
		$helperexpr
		eventrates = Array{Float64}(undef,nevents)
		$eventexprblock
		eventrates
	end
	end )

	# eval( quote
	# function fneventrates!(eventrates, t, A, interpdict, nevents) 
	# 	eventrates = fneventrates(t,A,interpdict, nevents)
	# 	# $Aex 
	# 	# $assexpr 
	# 	# $paex
	# 	# $helperexpr
	# 	# $eventexprblock
	# 	nothing
	# end
	# end )
	
	
	function eventrate(u,p,t) 
# println("entered eventrate")
		# fneventrates!(eventrates, t, A, interpdict) # A, interpdict 
		# Base.invokelatest( fneventrates!, eventrates, t, A, interpdict, nevents )
		# eventrates = fneventrates( eventrates, t, A, interpdict, nevents )
		eventrates = Base.invokelatest( fneventrates,  t, A, interpdict, nevents )
# if  t < 1.5
# println( eventrates )
# # any(isnan.(eventrates)) &&  @bp
# println(A)
# println(t)
# println("*************")
# end
		eventrates[isnan.(eventrates)] .= 0.0 
		max(0., sum(eventrates ))
	end
	leventrate(u,p,t) = 0.0 #sum(eventrates)/10.0 # TODO make this a parameter 
	# heventrate(u,p,t) = 100000.0# sum(eventrates)*10.0
	function heventrate(u,p,t)
		# ers =  fneventrates( t,All2,interpdict, nevents) 
		# rate with all samples extant
		ers = Base.invokelatest( fneventrates, t, All2, interpdict, nevents )
		10*sum(ers) 
	end
	#rint(u,p,t) = mindeltat  # mindeltat based on model solution delta t # mst/1e5
	rint(u,p,t) = mediandeltat
# rint(u,p,t) =  min( max(mst/1e4, (ushs[ ixsampleheight ]-t)/10.0) , 1e-4/sum(eventrates))
	# rint(u,p,t) = (0==sum(eventrates)) ? mst/10000.0 : min(mst/n , 10.0/sum(eventrates))

	function aff_event!(int)
		# fneventrates!(eventrates, int.t, A, interpdict) # A, interpdict 
		# Base.invokelatest( fneventrates!, eventrates, int.t, A, interpdict )
		eventrate( int.u, (), int.t )
# println(" &&&&&&&&&&&&&&&&&&&&& entered affevent &&&&&&&&&&&&&&&&&")
# println(eventrates)
# println(int.t)
@bp
		we = sample( Weights(eventrates) )
		rxntype = eventtypes[ we ]
		rxn = eventrxns[we]
		if rxntype == COALESCENT
			A[ rxn.recipient ] -= 1
			push!(events, Event(COALESCENT
			  , int.t 
			  , rxn.source
			  , rxn.recipient )
			)
		elseif rxntype == MIGRATION
			A[ rxn.recipient ] -= 1
			A[ rxn.source ] += 1
			push!(events, Event(MIGRATION
			  , int.t 
			  , rxn.source
			  , rxn.recipient )
			)
		end
		# int.u[we] += 1
		int.u[1] = 0.0 
# println( rxntype )
# println( int.u )
# println( eventrate( int.u, (), int.t ) )
# println( leventrate( int.u, (), int.t ) )
# println( heventrate( int.u, (), int.t ) )
# println( rint( int.u, (), int.t ) )

		cou = rand()
	end




	function sampcondition(u,t,integrator )::Bool
		# t>0 && t in sampleheights
		t == ushs[ ixsampleheight ] 
	end
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
		eventrate( integrator.u, (), integrator.t)
# println("entered sampeffect")
# println(integrator.t)
# println(eventrates)
		ixsampleheight = min( length( ushs ), ixsampleheight + 1 ) 
	end
	sampcb = DiscreteCallback( sampcondition, sampaffect! )
	
# 	vrj = VariableRateJump( eventrate, aff_event!; 
# 		 lrate = leventrate 
# 		, urate  = heventrate
# 		, rateinterval = rint
# 		)
# 	dcoprob = DiscreteProblem( zeros(length(eventrates))
# 	 ,  (0.0, mst)
# 	 , (;)
# 	)
# 	jmps = [vrj] 
# 	coprob = JumpProblem( dcoprob, Coevolve(), jmps...
# 	; dep_graph = [1])
# @bp
# mindeltat
# Base.invokelatest( fneventrates, 0.0, All2, interpdict, nevents )
# 	s = solve( coprob, SSAStepper() , callback=sampcb, tstops=ushs ) # TODO make stepper adjustable
#
# 		eventcb = ContinuousCallback(eventcondition, eventaffect! )
	
	function coodes!(du, u, p, t)
		# fneventrates!(eventrates, t, A, interpdict) # A, interpdict 
		# Base.invokelatest( fneventrates!, eventrates, t, A, interpdict )
		eventrate(u,p,t)
		eventrates[isnan.(eventrates)] .= 0.0 
		du[1] = max(0., sum(eventrates) ); 
	end


	function eventcondition(u,t,integrator)::Real
		#u[1] - cou
		exp(-u[1]) - cou # note: not integrator.u[1]
	end

	eventcb = ContinuousCallback(eventcondition, aff_event! )
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
	   , callback = cbs, tstops = ushs)


@bp 
	SimTree( events, model )
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
	
	assexprs = [ :( $(Symbol(v)) = interpdict[$v](t)) for (i,v) in enumerate( vcat(model.demes, model.nondemes) ) ]
	assexpr = Expr( :block, assexprs... )
	
	mparameters = model.parameters
	mparameters["mst"] = mst # NOTE copying here so the variable is accessible inside ODE expressions 
	paex =  Expr( :block, [ ( :($(Symbol(k)) = $v) ) for (k,v) in mparameters]... )
	helperexpr = Expr( :block, model.helperexprs... )
	Aex = Expr( :block, [ ( :($(Symbol("A_"*k)) = A[$k]) ) for k in keys(A)]... )

	function _sanitize_time_variable( expr )
		#= 
		if time (t) appears in expr, replace with (mst-t) for coalescent equations
				example input:
				:(if t > 80.6
					mu * hostB3
				else
					0.0
				end)
		=# 
		MacroTools.postwalk( expr ) do x 
			x == :t ? :(mst-t) : x 
		end
	end
	function _corateex( rxn )
		@assert rxn.type == RXN_BIRTH
		@assert !isnothing( rxn.recipient )
		src = rxn.source 
		rec = rxn.recipient

		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)
		Yₛ = :( max( $Aₛ,   $(Symbol(src))))
		Yᵣ = :( max( $Aᵣ,   $(Symbol(rec))))

		if src == rec 
			aex = :( $Aₛ * max(0.0,$Aᵣ-1.0) / ($Yₛ * $Yᵣ) )
		else
			aex = :( $Aₛ * $Aᵣ / ($Yₛ * $Yᵣ) )
		end
		Expr( :call, :*, aex,  _sanitize_time_variable( rxn.expr ) )
	end
	function _birthmigrateex( rxn )
		# forward time transmission s -> r, reverse time migration r -> s
		@assert rxn.type == RXN_BIRTH
		@assert !isnothing( rxn.recipient )
		src = rxn.source 
		rec = rxn.recipient

		Aₛ = Symbol("A_"*src)
		Aᵣ = Symbol("A_"*rec)
		Yₛ = :( max( $Aₛ,   $(Symbol(src))))
		Yᵣ = :( max( $Aᵣ,   $(Symbol(rec))))

		pex = :( clamp( ($Yₛ - $Aₛ) / $Yₛ , 0., 1. ) )
		aex = Expr( :call, :*, :( $Aᵣ / $Yᵣ ), pex )

		Expr( :call, :*, aex,  _sanitize_time_variable( rxn.expr )  )
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

	# function fneventrates!(eventrates, t, A, interpdict) 
	# 	eval( quote
	# 		$Aex 
	# 		$assexpr 
	# 		$paex
	# 		$helperexpr
	# 		$eventexprblock
	# 		nothing
	# 	end )
	# end 
	

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
		# fneventrates!(eventrates, t, A, interpdict) # A, interpdict 
		Base.invokelatest( fneventrates!, eventrates, t, A, interpdict )
		eventrates[isnan.(eventrates)] .= 0.0 
# println("~~~")
# println( A)
# println( [t, eventrates...] )
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
		# fneventrates!(eventrates, integrator.t, A, interpdict) # A, interpdict 
		Base.invokelatest( fneventrates!, eventrates, integrator.t, A, interpdict )
		eventrates[isnan.(eventrates)] .= 0.0 
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
	
# @bp
# DataFrame( events )|>print

	SimTree( events, model )
end

