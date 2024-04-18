module Coalescent

export SimTree, toRphylo, tonewick
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

const RXN_BIRTH = 0 
const RXN_MIG = 1 
const RXN_DEATH = 2 
const RXN_DYNVAR = 3 

const SAMPLE = 0
const COALESCENT = 1
const MIGRATION = 2
const RECOMBINATION = 3

const ALGO_MARKOV = "markov"
const ALGO_STATIONARY = "stationary" 


include("m0.jl")
include("a0.jl")
include("p0.jl")
include("s0.jl")

end;


using .Coalescent
using YAML 
using Test
using Plots
using Debugger 

conf = YAML.load_file("../test/q1.yaml")
m = ModelFGY( "../test/q1.yaml" )
o = solveodes( m )
plot( o )

s = SampleConfiguration( "../test/q1.yaml" )

# @run t = SimTree( m, s )
t = SimTree( m, s )
