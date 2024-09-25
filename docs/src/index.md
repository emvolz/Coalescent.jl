# Coalescent.jl

## Index

```@index
```

```@docs
Coalescent.Coalescent
```

## Functions

```@docs
SimTree
tonewick
ModelFGY
SampleConfiguration
solveodes
```

## Detailed Function Documentation

### SimTree

```@docs
SimTree(Ne::Float64, n::Int64)
SimTree(Ne::Float64, sampletimes::Array{Float64}, p...;  tmrcaguess::Union{Nothing,Float64}=nothing, algorithm=ALGO_MARKOV)
SimTree(Ne::Function, sampletimes::Array{Float64}, p...; tmrcaguess::Union{Nothing,Float64}=nothing, algorithm=ALGO_STATIONARY)
SimTree(Ne::Function, n::Int64, p...; tmrcaguess::Union{Nothing,Float64}=nothing, algorithm=ALGO_MARKOV)
SimTree(model::ModelFGY, sample::SampleConfiguration; computedescendants = false)
```

### ModelFGY

```@docs
ModelFGY(; modelname::String, birthrxn::Array{Reaction}, migrationrxn::Array{Reaction}, deathrxn::Array{Reaction}, nondemerxn::Array{Reaction}, demes::Array{String}, nondemes::Union{Nothing,Array{String}}, initial::Dict{String,Float64}, t0::Float64, tfin::Float64, parameters::Dict{String,Float64}, helperexprs::Union{Nothing,Array{Expr}})
ModelFGY(conffn::String)
```

### SampleConfiguration

```@docs
SampleConfiguration(conffn::String)
SampleConfiguration(; confstr::String)
```

### Other Functions

```@docs
tonewick(o)
Coalescent.distancematrix(t)
solveodes(model::ModelFGY; odemethod = :(Rosenbrock23()), res::Union{Missing,Int64} = missing)
```
