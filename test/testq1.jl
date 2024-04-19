using Coalescent
using YAML 
using Test
using Plots

conf = YAML.load_file("q1.yaml")
m = ModelFGY( "q1.yaml" )
o = solveodes( m )
plot( o )

s = SampleConfiguration( "q1.yaml" )

t = SimTree( m, s )

#= 
# TODO bug early termination sometimes and giving multifurc
 Warning: Instability detected. Aborting
└ @ SciMLBase ~/.julia/packages/SciMLBase/Avnpi/src/integrator_interface.jl:611
Simulated coalescent tree with 24 tips and 19 internal nodes
=#
