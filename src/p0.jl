# using Coalescent
using Statistics
using Random
using Distributions
using OrdinaryDiffEq
using Interpolations
using Plots

# """
#     Event
#
# Represents an event in the coalescent process.
#
# # Fields
# - `type::Int`: Type of the event (SAMPLE, COALESCENT, MIGRATION, or RECOMBINATION)
# - `height::Float64`: Time of the event
# - `source::Union{Nothing,String}`: Source deme of the event (if applicable)
# - `sink::Union{Nothing,String}`: Sink deme of the event (if applicable)
# """
mutable struct Event
	type::Int
	height::Float64
	source::Union{Nothing,String}
	sink::Union{Nothing,String}
end
# """
#     Event(t::Int, b::Real)
#
# Constructor for Event with only type and height specified.
#
# # Arguments
# - `t::Int`: Type of the event
# - `b::Real`: Height (time) of the event
#
# # Returns
# - `Event`: The constructed event
# """
function Event(t::Int,b::Real) 
	@assert t in [SAMPLE,COALESCENT,MIGRATION,RECOMBINATION] 
	Event( t, b, nothing, nothing )
end



"""
    SimTree

Represents a simulated coalescent tree.

# Fields
- `parent::Array{Union{Nothing,Int}}`: Parent nodes
- `child::Array{Union{Nothing,Int}}`: Child nodes
- `n::Int`: Number of tips
- `nNode::Int`: Number of internal nodes
- `edgelength::Array{Union{Nothing,Float64}}`: Edge lengths
- `heights::Array{Union{Nothing,Float64}}`: Node heights
- `tiplabs::Array{Union{Nothing,String}}`: Tip labels
- `shs::Array{Union{Nothing,Float64}}`: Sample heights
- `descendants::Union{Nothing,BitMatrix}`: Descendant matrix
- `daughters::Union{Nothing,Vector{Tuple{Int, Int, Int, Int}}}`: Daughter nodes
- `demes::Union{Nothing,Vector{String}}`: Demes for each node
"""
@kwdef mutable struct SimTree
	parent::Array{Union{Nothing,Int}}
	child::Array{Union{Nothing,Int}}
	n::Int
	nNode::Int
	edgelength::Array{Union{Nothing,Float64}}
	heights::Array{Union{Nothing,Float64}} # n + nNode 
	tiplabs::Array{Union{Nothing,String}}
	shs::Array{Union{Nothing,Float64}} 
	descendants::Union{Nothing,BitMatrix} = nothing # nNode * n 
	daughters::Union{Nothing,Vector{Tuple{Int, Int, Int, Int}}} = nothing # nNode elements 
	demes::Union{Nothing,Vector{String}} = nothing
end
Base.show( io::IO, x::SimTree) = print("""
Simulated coalescent tree with $(x.n) tips and $(x.nNode) internal nodes
Tip labels:
	$(x.tiplabs[1:min(10,x.n)]) ...

Rooted; includes branch lengths
""")


"""
    SimTree(events::Array{Event})::SimTree

Construct a SimTree from an array of Events.

# Arguments
- `events::Array{Event}`: Array of coalescent events

# Returns
- `SimTree`: The constructed simulated tree
"""
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
	daughters = Vector{Tuple{Int, Int, Int, Int }}(undef, nNode  )
	
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
			iu,iv = sample( 1:length(extant), 2, replace=false ) # as fast as using rand()*length|>floor|>+1|>Int 
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
			daughters[a-n] = (u,v,ie-1,ie)
			ie += 1 
			a -= 1
		else
			Error( "Not Implemented" )
		end
	end
	
	SimTree( parent, child, n, nNode, edgelength, heights, tiplabs,  heights[1:n]
	, nothing, daughters # descendants , daughters 
	, nothing # deme
	)
end


# """
#     _sim_markov(Ne::Function, sampletimes::Array{Float64}, tmrcaguess::Float64, p...)::SimTree
#
# Simulate a coalescent tree using a Markovian coalescent algorithm.
#
# # Arguments
# - `Ne::Function`: Effective population size function
# - `sampletimes::Array{Float64}`: Array of sample times
# - `tmrcaguess::Float64`: Initial guess for the time to most recent common ancestor
# - `p...`: Additional parameters for the Ne function
#
# # Returns
# - `SimTree`: The simulated coalescent tree
# """
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

	function coodes!(du, u, pars, t)
		N = Ne(t, p... )
		du[1] = max(0., A*(A-1.)/(2N)) # cumulative hazard 
	end

	# TODO optimise samp callback by pretabulating sample time and sample sizes
	function sampcondition(u,t,integrator )::Bool
		t>0 && t in sampleheights
	end
	function sampaffect!(integrator)
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
	integ = AutoTsit5(Rosenbrock23())
	# integ = RadauIIA5() 
	# integ = Tsit5()

	pr = ODEProblem( coodes!
		 , [0.0]
		 #, [ 1.0 ]
		 , ( 0., tmrcaguess )
		)
	s = solve( pr, integ 
	   , callback = cbs, tstops = unique( sampleheights ) )
	u = s.u 
	taxis = s.t
	tfin0 = s.t[end] 
	dxt = s.t[end] - s.t[end-1]
	tfin = tfin0 
	while cocount < (n-1)
		pr = ODEProblem( coodes!
		  , s.u[end]
		  , ( s.t[end], tfin )
		)
		s = solve( pr, integ
		   , callback = cbs
		   , tstops = unique( sampleheights ) 
		   # , dt = dxt 
		   # , dtmin = dxt
		   # , force_dtmin = true 
		   )
		u = [u ; s.u ]
		taxis = [taxis ; s.t ]
		tfin += tfin0

		# N = Ne(s.t[end], p... )
		# @show last(u, 5 )
		# @show last( taxis, 5) 
		# @show cocount
		# @show tfin 
		# @show A 
		# @show N
		# @show max(0., A*(A-1.)/(2Ne(s.t[end], p... ))) 
		# @show dxt 
		# @assert !isnan(N)
	end
	
	tr = SimTree( [[Event(COALESCENT, h) for h in coheights] ; 
		[Event(SAMPLE, h) for h in sampleheights]
	])
	tr

end



"""
    SimTree(Ne::Function, sampletimes::Array{Float64}, p...; tmrcaguess::Union{Nothing,Float64}=nothing, algorithm=ALGO_STATIONARY)::SimTree

Simulate a coalescent tree with flexible Ne function and sampling times.

# Arguments
- `Ne::Function`: Effective population size function
- `sampletimes::Array{Float64}`: Array of sample times
- `tmrcaguess::Float64`: Initial guess for the time to most recent common ancestor
- `p...`: Additional parameters for the Ne function

# Keywords
- `algorithm::String = ALGO_MARKOV`: Algorithm to use for simulation (ALGO_STATIONARY or ALGO_MARKOV)

# Returns
- `SimTree`: The simulated coalescent tree
"""
function SimTree(Ne::Function, sampletimes::Array{Float64}, p...
	;tmrcaguess::Union{Nothing,Float64}=nothing,  algorithm=ALGO_MARKOV)::SimTree
	@assert algorithm in [ALGO_STATIONARY, ALGO_MARKOV]
	tmrcaguess = isnothing(tmrcaguess) ? Ne(0.0,p...) : tmrcaguess
	if algorithm == ALGO_MARKOV
		return _sim_markov(Ne, sampletimes, tmrcaguess, p...)
	end
	if algorithm == ALGO_STATIONARY 
		return _sim_stationary( Ne, sampletimes, tmrcaguess, p... )
	end
end

function _sim_stationary( Ne::Function, sampletimes::Array{Float64}, tmrcaguess::Float64 , p... )
	print( """Simulating coalescent, sample size =$(length(sampletimes))
Initial guess of time to most recent common ancestor (before most recent sample): $(tmrcaguess)
Deterministic approximation to coalescent time distribution
User-specified Ne(t) function 
""")
	@warn "Entering development code"
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
		N = Ne(t, p...) 
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
		sampadds += 1
		integrator.u[1] += sum( integrator.t .== sampleheights )
	end
	sampcb = DiscreteCallback( sampcondition, sampaffect! )
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
	@show u

	cointerp = linear_interpolation( sort!([ uu[2] for uu in u ] ./ float(n-1)) , taxis )
	g() = sort!( cointerp.( rand(n-1))  )

	# ensure cotimes compatible with samp times 
	cou = sort!(rand(n-1))
	try
		X = [cou cointerp.(cou) -ones(n-1) repeat([COALESCENT], n-1) ]
	catch
		Error("Interpolation Error.")
	end
	X = [X ; 
	    repeat([missing],n) sampleheights ones(n) repeat([SAMPLE], n) ]
	X = X[ sortperm( X[:,2] ), :]
	A = cumsum( X[:,3] )
	while any(A .< 1.)
		i = findfirst( A .< 1. )
		_cou = X[i,1] + (1.0-X[i,1])*rand()
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
#= 
"""
    SimTree(Ne::Float64, n::Int64, tmrcaguess::Float64, p...; algorithm=ALGO_STATIONARY)::SimTree

Simulate a coalescent tree with constant Ne and n samples at time 0.

# Arguments
- `Ne::Float64`: Constant effective population size
- `n::Int64`: Number of samples
- `tmrcaguess::Float64`: Initial guess for the time to most recent common ancestor
- `p...`: Additional parameters (unused for constant Ne)

# Keywords
- `algorithm::String = ALGO_STATIONARY`: Algorithm to use for simulation (ALGO_STATIONARY or ALGO_MARKOV)

# Returns
- `SimTree`: The simulated coalescent tree
"""
function SimTree(Ne::Float64, n::Int64, tmrcaguess::Float64, p... ; algorithm=ALGO_STATIONARY)::SimTree
	_Ne(t,p...) = Ne
	SimTree( _Ne 
	 , repeat( [0.], n) 
	 , tmrcaguess
	 , p...
	 ; algorithm
	)
end =#

"""
    SimTree(Ne::Float64, sampletimes::Array{Float64}, p...;  tmrcaguess::Union{Nothing,Float64}=nothing, algorithm=ALGO_MARKOV)::SimTree

Simulate a coalescent tree with constant Ne and flexible sampling times.

# Arguments
- `Ne::Float64`: Constant effective population size
- `sampletimes::Array{Float64}`: Array of sample times
- `tmrcaguess::Float64`: Initial guess for the time to most recent common ancestor
- `p...`: Additional parameters (unused for constant Ne)

# Keywords
- `algorithm::String = ALGO_STATIONARY`: Algorithm to use for simulation (ALGO_STATIONARY or ALGO_MARKOV)

# Returns
- `SimTree`: The simulated coalescent tree
"""
function SimTree(Ne::Float64, sampletimes::Array{Float64}, p... ;  tmrcaguess::Union{Nothing,Float64}=nothing, algorithm=ALGO_MARKOV)::SimTree
	_Ne(t,p...) = Ne
	SimTree( _Ne
	 , sampletimes 
	 , p...
 	 ; tmrcaguess,  algorithm
	)
end

"""
    SimTree(Ne::Function, n::Int64,p...; tmrcaguess::Union{Nothing,Float64}=nothing, algorithm=ALGO_MARKOV)::SimTree

Simulate a coalescent tree with flexible Ne function and n samples at time 0.

# Arguments
- `Ne::Function`: Effective population size over time function
- `n::Int64`: Number of samples
- `tmrcaguess::Float64`: Initial guess for the time to most recent common ancestor
- `p...`: Additional parameters for the Ne function

# Keywords
- `algorithm::String = ALGO_MARKOV`: Algorithm to use for simulation (ALGO_STATIONARY or ALGO_MARKOV)

# Returns
- `SimTree`: The simulated coalescent tree
"""
function SimTree(Ne::Function, n::Int64,  p... ; tmrcaguess::Union{Nothing,Float64}=nothing, algorithm=ALGO_MARKOV)::SimTree
	tmrcaguess = isnothing(tmrcaguess) ? Ne(0.0,p...) : tmrcaguess 	
	SimTree( Ne
	 , repeat( [0.], n) 
	 , p...
 	 ; tmrcaguess
 	 , algorithm
	)
end

"""
    SimTree(Ne::Float64, n::Int64)::SimTree

Simulate a coalescent tree with constant Ne and n samples at time 0, without specifying tmrcaguess.

# Arguments
- `Ne::Float64`: Constant effective population size
- `n::Int64`: Number of samples

# Returns
- `SimTree`: The simulated coalescent tree
"""
function SimTree( Ne::Float64, n::Int64)::SimTree 
	@assert n > 1
	@assert Ne > 0 
	events = Array{Event}(undef, n+n-1)
	for i in 1:n
		events[i] = Event(SAMPLE, 0.0 )
	end
	inis = Array{Float64}(undef, n-1 )
	for i in n:-1:2
		inis[n-i+1] = rand( Exponential( 2Ne/ (i*(i-1))))
	end
	nhs = cumsum( inis )
	for (i, h) in enumerate( nhs )
		events[n+i] = Event( COALESCENT, h )
	end 
	SimTree( events )
end 

#= function toRphylo(stre)
	edge = [stre.parent stre.child]
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
=#



# TODO compute mdt, dgtrs in SimTree 
"""
    tonewick(o)

Convert a SimTree to a Newick format string.

# Arguments
- `o::SimTree`: The SimTree to convert

# Returns
- `String`: Newick format representation of the tree
"""
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

"""
    distancematrix(t)::Matrix{Float64}

Compute the distance matrix for a SimTree.

# Arguments
- `t::SimTree`: The SimTree to compute distances for

# Returns
- `Matrix{Float64}`: The computed distance matrix
"""
function distancematrix(t)::Matrix{Float64}
	(isnothing(t.descendants)) && throw("SimTree must be generated with computedescendants=true.")

	d = zeros( t.n, t.n )
	dnt = zeros( t.n + t.nNode, t.n )

	# TODO should not rely on heights for postorder 
	ponodes = sortperm( t.heights[t.n+1:end] ) .+ t.n
	for a in ponodes
		u,v,ieu,iev = t.daughters[a-t.n] 
		elu = t.edgelength[ ieu ]
		elv = t.edgelength[ iev ]
		# descu = u>t.n ? findall( t.descendants[u-t.n,:] ) : u 
		# descv = v>t.n ? findall( t.descendants[v-t.n,:] ) : v 

		# define desc tips of dautghter and fill in distance from a to these desc
		if u > t.n 
			descu = findall( t.descendants[u-t.n,:] )
			dnt[ a, descu ] .= elu .+ dnt[u, descu ]
		else 
			descu = u 
			dnt[ a, descu ] = elu + dnt[u, descu ]
		end
		
		# define desc tips of dautghter and fill in distance from a to these desc
		if v > t.n 
			descv = findall( t.descendants[v-t.n,:] )
			dnt[ a, descv ] .= elv .+ dnt[v, descv ]
		else 
			descv = v
			dnt[ a, descv ] = elv + dnt[v, descv ]
		end

		# fill in descu * descv elements of distance mat 
		if (u > t.n) & (v > t.n)
			d[ descu, descv ] .= dnt[ a, descv ]' .+ dnt[ a, descu ]
			d[ descv, descu ] .= d[ descu, descv ]'
		elseif  (u > t.n) & (v <= t.n) 
			d[ descu, descv ] .= dnt[ a, descu ] .+ dnt[ a, descv ]
			d[ descv, descu ] .= d[ descu, descv ]
		elseif (u <= t.n) & (v > t.n) 
			d[ descu, descv ] .= dnt[ a, descu ] .+ dnt[ a, descv ]
			d[ descv, descu ] .= d[ descu, descv ]
		else 
			d[ descu, descv ] = dnt[ a, descv ] + dnt[ a, descu ]
			d[ descv, descu ] = d[ descu, descv ]
		end

		# for uu in t.descendants[u] 
		# 	d[uu, t.descendants[v]] .= dnt[a, t.descendants[v]] .+ dnt[a, uu] 
		# end
	end

	d
end
