using Coalescent
using Test

function _Ne( t, p... )
	# println( "_Ne $(p)" )
	r, N0 = p[1] 
	exp( -r*t ) * N0 
end

function _Ne2( t ) # dnw 
	# println( "_Ne $(p)" )
	r = 1.
	N0 = 1e3
	exp( -r*t ) * N0 
end


@testset "Coalescent.jl" begin
    @test begin
    	jtr = SimTree( 1.0, Int(1e4), 10.,  1., 1e3, algorithm="markov" ) 
    	print( jtr )
    	true
    end
end


