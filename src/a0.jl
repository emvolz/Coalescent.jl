#= 
#sample configuration
=# 

using YAML

struct SampleConfiguration 
	deme::Array{String}
	time::Array{Float64}
end 
function SampleConfiguration( conffn::String )
	conf = YAML.load_file( conffn )
	samps = conf["sample"] 

	deme = Array{String}()
	time = Array{Float}()

	for s in samps 
		if "table" ∈ keys(s)
			...
		elseif "deme" ∈ keys(s)
			if typeof( s["deme"] ) == String 
				push!( deme, s["deme"] )
			elseif typeof(s["deme"] ) == Array 
			end
		end
	end 
end
