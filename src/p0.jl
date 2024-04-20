# using Coalescent
using Statistics
using Random
using Distributions
using RCall
using OrdinaryDiffEq
using Interpolations
using Plots
using Debugger

struct Event
	type::Int
	height::Float64
	source::Union{Nothing,String}
	sink::Union{Nothing,String}
end
function Event(t::Int,b::Real) 
	@assert t in [SAMPLE,COALESCENT,MIGRATION,RECOMBINATION] 
	Event( t, b, nothing, nothing )
end



mutable struct SimTree
	parent::Array{Union{Nothing,Int}}
	child::Array{Union{Nothing,Int}}
	n::Int
	nNode::Int
	edgelength::Array{Union{Nothing,Float64}}
	heights::Array{Union{Nothing,Float64}}
	tiplabs::Array{Union{Nothing,String}}
end
Base.show( io::IO, x::SimTree) = print("""
Simulated coalescent tree with $(x.n) tips and $(x.nNode) internal nodes
Tip labels:
	$(x.tiplabs[1:min(10,x.n)]) ...

Rooted; includes branch lengths
""")
function SimTree(events::Array{Event})::SimTree
	ix = sortperm( [ e.height for e in events ] )
	events = events[ ix ]
	# ses = [ e for e in events if e.type == SAMPLE ]
	# n = length( ses )
	n = sum( [ e.type == SAMPLE for e in events ] )
	# coes = [ e for e in events if e.type == COALESCENT ]
	nNode = sum( [e.type == COALESCENT for e in events ] )
	# nNode = length( coes )

	tiplabs = [ "t.$(i)" for i in 1:n ]
	nodes = range( 1, n + nNode ) |> collect 
	
	nedges = n + nNode - 1
	parent = Array{Union{Nothing,Int}}(nothing, nedges )
	child = Array{Union{Nothing,Int}}(nothing, nedges )
	edgelength = Array{Union{Nothing,Float64}}(nothing, nedges )
	heights = Array{Union{Nothing,Float64}}( nothing, n + nNode )
	
	A=0
	ie = 1
	extant = Array{Int}(undef, 0)
	sizehint!( extant, n + nNode )
	a = n + nNode # internal node counter, goes down 
	su = 1 # tip counter, goes up 
	
	iu = 0 
	iv = 0

	for e in events
		@assert isnothing( e.source ) # unstructured models only 
		if e.type == SAMPLE
			A+=1 
			heights[su] = e.height 
			push!( extant, su )
			su += 1
		elseif e.type == COALESCENT
			heights[a] = e.height
			try
				iu,iv = sample( 1:length(extant), 2, replace=false )
			catch
				@bp
			end
			u,v = extant[ [iu,iv] ]
			deleteat!( extant, sort([iu,iv]) )
			push!(extant, a )
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
		else
			Error( "Not Implemented" )
		end
	end
	
	SimTree( parent, child, n, nNode, edgelength, heights, tiplabs )
end

function _sim_markov(Ne::Function
		     , sampletimes::Array{Float64}
		     , tmrcaguess::Float64
		     , p... 
	)::SimTree
	print( """Simulating coalescent, sample size =$(length(sampletimes))
Initial guess of time to most recent common ancestor (before most recent sample): $(tmrcaguess)
Markovian coalescent algorithm
User-specified Ne(t) function 
""")
	n = length( sampletimes )
	@assert n > 1
	@assert tmrcaguess > 0. 
	@assert diff( collect( extrema( sampletimes)))[1] < tmrcaguess 
	sort!(sampletimes)
	mst = maximum( sampletimes )
	sampleheights = reverse(  mst .- sampletimes ) 
	pars = Tuple( collect(p) )
	
	sampadds = 0 
	cou = rand() 
	sampcount = 0 
	cocount = 0 
	coheights = Array{Float64}(undef,0); sizehint!(coheights, n-1)
	A = sum(sampleheights .== 0.) |> float 

	function coodes!(du, u, p, t)
		N = Ne(t, pars )
		# du[1] = -u[1] * A*(A-1.0)/(2N)
		du[1] = max(0., A*(A-1.)/(2N)) # cumulative hazard 
	end

	# TODO optimise samp callback by pretabulating sample time and sample sizes
	function sampcondition(u,t,integrator )::Bool
		t>0 && t in sampleheights
	end
	function sampaffect!(integrator)
		 #TODO does this work with multiple samples at integrator.t? probably not
		sampadds += 1
		# A += 1.
		A += sum( sampleheights .== integrator.t  )
	end
	sampcb = DiscreteCallback( sampcondition, sampaffect! )
	
	function cocondition(u,t,integrator)::Real
		#u[1] - cou
		exp(-u[1]) - cou # note: not integrator.u[1]
	end
	function coaffect!(integrator)
		cocount+=1 
		A-=1. 
		# integrator.u[1] = 1.
		integrator.u[1] = 0. 
		cou = rand()
		push!(coheights, integrator.t)
	end
	cocb = ContinuousCallback( cocondition, coaffect! )

	cbs = CallbackSet( sampcb, cocb )
	
	# integ = Rosenbrock23()
	# integ = AutoTsit5(Rosenbrock23())
	integ = Tsit5()

	pr = ODEProblem( coodes!
		 , [0.0]
		 #, [ 1.0 ]
		 , ( 0., tmrcaguess )
		)
	s = solve( pr, integ 
	   , callback = cbs, tstops = unique( sampleheights ) )
	u = s.u 
	taxis = s.t 
	tfin = s.t[end]*1.5
	while cocount < (n-1)
		pr = ODEProblem( coodes!
		  , s.u[end]
		  , ( s.t[end], tfin )
		)
		s = solve( pr, integ
		   , callback = cbs, tstops = unique( sampleheights ) )
		u = [u ; s.u ]
		taxis = [taxis ; s.t ]
		tfin = s.t[end]*1.5
	end
	
	tr = SimTree( [[Event(COALESCENT, h) for h in coheights] ; 
		[Event(SAMPLE, h) for h in sampleheights]
	])
	tr

end

function SimTree(Ne::Function, sampletimes::Array{Float64}, tmrcaguess::Float64
		 , p...
	; algorithm=ALGO_STATIONARY)::SimTree

	@assert algorithm in [ALGO_STATIONARY, ALGO_MARKOV]
	if algorithm == ALGO_MARKOV
		return _sim_markov(Ne, sampletimes, tmrcaguess, p...)
	end
	print( """Simulating coalescent, sample size =$(length(sampletimes))
Initial guess of time to most recent common ancestor (before most recent sample): $(tmrcaguess)
Deterministic approximation to coalescent time distribution
User-specified Ne(t) function 
""")
	n = length( sampletimes )
	@assert n > 1
	@assert tmrcaguess > 0. 
	@assert diff( collect( extrema( sampletimes)))[1] < tmrcaguess 
	sort!(sampletimes)
	mst = maximum( sampletimes )
	sampleheights = reverse(  mst .- sampletimes ) 
	pars = Tuple( collect(p) )
	# @show pars 
	function lttode!(du, u, p, t) 
		# println( "lttode $(p)" )
		N = Ne(t, pars) 
		du[1] = -u[1]*(u[1]-1.) / (2N)
		du[2] = -du[1] 
	end
	pr = ODEProblem( lttode!
		 , [ sum(sampleheights .== 0.), 0. ]
		 , ( 0., tmrcaguess )
		)
	function sampcondition(u,t,integrator )
		o = t>0 && t in sampleheights
		# println("integ time $(t)")
		# if o
		# 	println( "sampcond reached $(t)" )
		# end
		o
	end
	sampadds = 0 
	function sampaffect!(integrator)
		# println( "sampaffect $(integrator.u)" )
		# println( "sampaffect $(integrator.t)" )
		# println( "sampaffect $(integrator.t .== sampleheights)" )
		# println( "sampaffect $(sum(integrator.t .== sampleheights))" )
		sampadds += 1
		integrator.u[1] += sum( integrator.t .== sampleheights )
	end
	sampcb = DiscreteCallback( sampcondition, sampaffect! )
	# s = solve( pr, AutoTsit5(Rosenbrock23())
	#    , callback = sampcb, tstops = unique( sampleheights ) )
	s = solve( pr, Rosenbrock23()
	   , callback = sampcb, tstops = unique( sampleheights ) )
	u = s.u 
	taxis = s.t 
	afin = s.u[end][1]
	tfin = s.t[end]*1.5
	while afin > 1.1
		pr = ODEProblem( lttode!
		  , s.u[end]
		  , ( s.t[end], tfin )
		)
		s = solve( pr, AutoTsit5(Rosenbrock23())
		   , callback = sampcb, tstops = unique( sampleheights ) )
		u = [u ; s.u ]
		taxis = [taxis ; s.t ]
		afin = s.u[end][1]
		tfin = s.t[end]*1.5
	end
	@show sampadds 

	cointerp = linear_interpolation( sort!([ uu[2] for uu in u ] ./ float(n-1)) , taxis )
	g() = sort!( cointerp.( rand(n-1))  )

	# ensure cotimes compatible with samp times 
	cou = sort!(rand(n-1))
	try
		X = [cou cointerp.(cou) -ones(n-1) repeat([COALESCENT], n-1) ]
	catch
		@bp 
		Error("Interpolation Error.")
	end
	X = [X ; 
	    repeat([missing],n) sampleheights ones(n) repeat([SAMPLE], n) ]
	X = X[ sortperm( X[:,2] ), :]
	A = cumsum( X[:,3] )
	while any(A .< 1.)
		i = findfirst( A .< 1. )
		_cou = X[i,1] + (1.0-X[i,1])*rand()
		# @show "########"
		# @show A
		# @show i 
		# @show X[i,:]
		println("Resampling coalescent times, index $(i), height $(X[i,2])")

		X[i,1] = _cou 
		X[i,2] = cointerp(_cou)
		X = X[ sortperm( X[:,2] ), :]
		A = cumsum( X[:,3] )
	end
	
	# @bp 
	events = [ Event( Int(x), h ) for (x,h) in zip( X[:,4], X[:,2] ) ]
	
	tr = SimTree( events )
	tr
end

function SimTree(Ne::Float64, n::Int64, tmrcaguess::Float64, p... ; algorithm=ALGO_STATIONARY)::SimTree
	_Ne(t,p...) = Ne
	SimTree( _Ne 
	 , repeat( [0.], n) 
	 , tmrcaguess
	 , p...
	 ; algorithm
	)
end

function SimTree(Ne::Float64, sampletimes::Array{Float64}, tmrcaguess::Float64, p... ; algorithm=ALGO_STATIONARY)::SimTree
	_Ne(t,p...) = Ne
	SimTree( _Ne
	 , sampletimes 
	 , tmrcaguess
	 , p...
 	 ; algorithm
	)
end

function SimTree(Ne::Function, n::Int64, tmrcaguess::Float64, p... ; algorithm=ALGO_STATIONARY)::SimTree
	SimTree( Ne
	 , repeat( [0.], n) 
	 , tmrcaguess
	 , p...
 	 ; algorithm
	)
end

function toRphylo(stre)
	edge = [stre.parent stre.child]
	# @rput edge stre.edgelength stre.n stre.nNode stre.tiplabs
	# R"""
	# library( ape )
	# edge <- unlist( edge )|> matrix( ncol =2 )
	# tr <- structure(list(n = n, Nnode = nnodes, edge = edge, edge.length=unlist(edgelength), tip.label=unlist(tiplabs)), class='phylo')
	# tr <- reorder(tr)
	# trnwk <- write.tree(tr)
	# 1
	# """;
	# @rget trnwk;
	R"""
	library(ape)
	tr <- structure(list(n = $(stre.n), Nnode = $(stre.nNode)
	  , edge = $(edge)
	  , edge.length=unlist( $(stre.edgelength) )
	  , tip.label=unlist($(stre.tiplabs))
	), class='phylo')
	tr <- ape::reorder.phylo(tr)
	trnwk <- ape::write.tree(tr)
	print( tr )
	1
	""";
	# @rget tr
	R"tr"
end




# TODO compute mdt, dgtrs in SimTree 
function tonewick(o)
	edge = [o.parent o.child]
	n = o.n 
	nnode = o.nNode 
	nwks = fill("", n + nnode)
	nwks[1:n] = [string(o.tiplabs[u], ":", o.edgelength[edge[:, 2] .== u][1]) for u in 1:n]
	queue = collect(1:n)

	# Compute max distance to tip
	mdt = fill(-Inf, n + nnode)
	mdt[1:n] .= 0
	edge0 = edge[sortperm(edge[:, 2]), :]
	anc = Array{Union{Nothing,Int}}(nothing, n+nnode) 
	anc[edge0[:, 2]] = edge0[:, 1]

	while !isempty(queue)
		for u in queue
			a = anc[u]
			if !isnothing(a)
				mdt[(a)] = max(mdt[(a)], mdt[u] + 1)
			end
		end
		queue = filter(!isnothing, anc[queue])
	end

	ord = setdiff(sortperm(mdt), 1:n)
	for a in ord
		dgtrs = edge0[edge0[:, 1] .== a, 2]
		nwks[a] = join(nwks[dgtrs], ",")
		edle = a in edge[:, 2] ? o.edgelength[edge[:, 2] .== a][1] : 0.0
		nwks[a] = "($(nwks[a])):$edle"
	end
	nwk = nwks[ord[end]] * ";"
	nwk
end


