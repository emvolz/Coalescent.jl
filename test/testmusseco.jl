using YAML 
using Test
using Plots
using OrdinaryDiffEq

conf = YAML.load_file( "test/musseco.yaml")

using Revise
using Coalescent

using Debugger 

#= @run SampleConfiguration( "test/musseco.yaml" ) 
 @run  ModelFGY( "test/musseco.yaml" ) 
=#
s = SampleConfiguration( "test/musseco.yaml" ) 
m = ModelFGY( "test/musseco.yaml" ) 

 # @run o =  solveodes( m )
o =  solveodes( m )
plot(o)

t = SimTree( m, s )
write("test.nwk" , tonewick(t) )
