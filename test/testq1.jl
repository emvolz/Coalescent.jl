# using Coalescent
using YAML 
using Test
using Plots

using Revise
using Coalescent

# include("../src/Coalescent.jl")
# using .Coalescent
# include("../src/s0.jl")

conf = YAML.load_file("test/q1.yaml")
m = ModelFGY( "test/q1.yaml" )
o = solveodes( m )
plot( o )

s = SampleConfiguration( "test/q1.yaml" )

using Debugger
# @run t = SimTree( m, s )
t = SimTree( m, s )

#= 
# TODO bug early termination sometimes and giving multifurc
 Warning: Instability detected. Aborting
└ @ SciMLBase ~/.julia/packages/SciMLBase/Avnpi/src/integrator_interface.jl:611
Simulated coalescent tree with 24 tips and 19 internal nodes
=#
