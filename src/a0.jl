#= 
#sample configuration
=# 

using YAML
using CSV
using DataFrames
using Distributions
using Random
# using Debugger 

struct SampleConfiguration 
	sconf::Array{Tuple{Union{Nothing,String}, Float64}}
end 
function SampleConfiguration( confstr::String )
	# conf = YAML.load_file( conffn )
	conf = YAML.load( confstr )
	samps = conf["sample"] 

	deme = Array{String}(undef,0)
	time = Array{Float64}(undef,0)

	for s in samps 
		if "table" ∈ keys(s)
			stx = DataFrame( CSV.File( s["table"] ))
			push!( time, stx.sample_time... )
			demes = [ String(x) for x in stx.deme ] # may need to promote from String1
			push!( deme, demes... )
		elseif "time" ∈ keys(s)
			if s["time"] isa String
				tx = ( Meta.parse( s["time"] ) ) |> eval |> collect 
			elseif s["time"] isa Number
				tx = s["time"] 
			else 
				throw( s )
			end
			ltx = length( tx )
			if ltx == 1 && "size" ∈ keys(s)
				@assert s["size"] isa Number 
				ltx = s["size"] 
				tx = repeat( [tx], ltx )
			end
			push!( time, tx... )
			if "deme" ∈ keys(s) 
				push!( deme, repeat( [s["deme"]], ltx )... )
			else 
				push!( deme, repeat( ["NULL"], ltx )... )
			end
		end
	end 
	SampleConfiguration( collect( zip( deme, time )) )
end



