# using Coalescent
using YAML 
using Test
using Plots
using OrdinaryDiffEq

using Revise
using Coalescent

fn = "test/sirexponential.yaml"
conf = YAML.load_file(fn)
m = ModelFGY( fn )
o = solveodes( m ; res = 100 )
# o1 = solveodes( m, integrator = RadauIIA5 ) 
# o2 = solveodes( m, integrator = Rodas4P) 
# plot( o )

s = SampleConfiguration( fn )

using Debugger
# @run t = SimTree( m, s )
t = SimTree( m, s )
