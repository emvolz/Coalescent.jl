using Test

using Revise
using Coalescent

using Debugger

@time t = SimTree( 1e4, 100_000 )
# @time t = SimTree( 1e4, 1_000_000 )
t isa SimTree 
println(t) 
