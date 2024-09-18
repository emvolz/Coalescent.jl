"""
Coalescent

Coalescent simulation and analysis with flexible markup of demographic processes and sampling patterns. 

This module provides tools for:

- Simulating coalescent trees
- Handling structured population models
- Solving ordinary differential equations for population dynamics
- Managing sample configuration

"""
module Coalescent

export SimTree,  tonewick
export ModelFGY, solveodes
export SampleConfiguration

using YAML 
using CSV
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

"Reaction types defining rates of genealogical events." 
const RXN_BIRTH = 0 
const RXN_MIG = 1 
const RXN_DEATH = 2 
const RXN_DYNVAR = 3 

"Events in a coalescent process." 
const SAMPLE = 0
const COALESCENT = 1
const MIGRATION = 2
const RECOMBINATION = 3

"Alternative models of the coalescent distribution. "
const ALGO_MARKOV = "markov"
const ALGO_STATIONARY = "stationary" 


include("m0.jl")
include("a0.jl")
include("p0.jl")
include("s0.jl")

end;
