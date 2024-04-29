# using Coalescent
using YAML 
using Test
using Plots
using OrdinaryDiffEq

using Revise
using Coalescent

conf = YAML.load_file("test/q1.yaml")
m = ModelFGY( "test/q1.yaml" )
o = solveodes( m ; res = 10000 )
# o1 = solveodes( m, integrator = RadauIIA5 ) 
# o2 = solveodes( m, integrator = Rodas4P) 
# plot( o )

s = SampleConfiguration( "test/q1.yaml" )

using Debugger
# @run t = SimTree( m, s )
t = SimTree( m, s )
