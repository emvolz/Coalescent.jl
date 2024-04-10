using Coalescent
using YAML 
using Test
using Plots

conf = YAML.load_file("q1.yaml")
m = ModelFGY( "q1.yaml" )
o = solveodes( m )
plot( o )


