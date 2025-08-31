"""
Coalescent

Coalescent simulation and analysis with flexible markup of demographic processes and sampling patterns. 

This module provides tools for:

- Simulating coalescent trees
- Handling structured population models  
- Solving ordinary differential equations for population dynamics
- Managing sample configuration
- DSL for defining models in pure Julia (NEW!)

"""
module Coalescent

export SimTree,  tonewick
export ModelFGY, solveodes
export SampleConfiguration
# DSL exports - the new treasure for model definition!
export @model, CoalescentModel, to_modelfgy

using YAML 
using CSV
using DataFrames
using OrdinaryDiffEq
using Statistics
using Random
using Distributions
using Interpolations
using Plots
using StatsBase

# "Reaction types defining rates of genealogical events." 
const RXN_BIRTH = 0 
const RXN_MIG = 1 
const RXN_DEATH = 2 
const RXN_DYNVAR = 3 

# "Events in a coalescent process." 
const SAMPLE = 0
const COALESCENT = 1
const MIGRATION = 2
const RECOMBINATION = 3

# "Alternative models of the coalescent distribution. "
const ALGO_MARKOV = "markov"
const ALGO_STATIONARY = "stationary" 


include("m0.jl")
include("a0.jl")
include("p0.jl")
include("s0.jl")
# Include our new DSL components - the crown jewels!
include("dsl.jl")
include("dsl_conversion.jl")

end;
