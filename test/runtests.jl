using Coalescent
using YAML 
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

yamltest_model = """
modelname: SEIR
version: 0.0.1
births:
  - source: I
    recipient: E
    rate: beta*S*I/N

migrations: 
  - source: E 
    recipient: I
    rate: gamma1*E

deaths:
  - deme: I
    rate: gamma2*I

parameters:
  - name: beta
    value: 3.0

  - name: gamma1
    value: 2.0

  - name: gamma2
    value: 2.0

# All variables which change through time. This should include all demes. 
dynamic_variables:
  - name: I
    initial_value: 1.0
  
  - name: E
    initial_value: 0.0 
  
  # Some dynamic variables are not demes (cannot be sampled). These should includ an `ode` to specify the rate of change over time
  - name: S
    initial_value: 1e5
    ode: -beta*S*I/N

  - name: R
    initial_value: 0.0
    ode: gamma2 * I

# Other variables which can be derived from other dynamic variables or parameters. These are not necessary, but can help simplify expressions
helpers: 
  - name: N
    definition: S + E + I + R

# time period over which simulation is carried out
time: 
  initial: 1.0
  final: 35.0 

# The number of simulations to carry out 
simulations: 10 

"""


yamltest_sampling = """
# Various ways to specify when and from which demes samples are collected
sample:
  # Take a single sample from deme I at time 10
  - deme: I 
    time: 10.
  # Take 10 samples from deme E at time 8
  - deme: E
    time: 8.0
    size: 10
  # Executable Julia code to define sample times. Takes 5 samples from deme I at equal intervals between times 5 and 10. 
  - deme: I
    time: range( 5., 10., 5 )
  # Another example of parsing and executing Julia code; sampling 10 times from a Normal distribution
  - deme: I
    time: rand( Normal(7.5), 10 )
  # Read sample time information from a table
  - table: ../inst/sampletimes.csv # A table with columns <sample_time>, <deme>

"""

@testset "Coalescent.jl" begin
	@test begin
		jtr = SimTree( 1.0, Int(1e4), 10.,  1., 1e3, algorithm="markov" ) 
		print( jtr )
		true
	end
	@test begin 
		SampleConfiguration( yamltest_sampling)
		true 
	end
end


