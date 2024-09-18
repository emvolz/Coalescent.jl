#= 
Sample configuration
=# 

using YAML
using CSV
using DataFrames
using Distributions
using Random

"""
    SampleConfiguration

Represents a configuration for sampling in coalescent simulations.

# Fields
- `sconf::Array{Tuple{Union{Nothing,String}, Float64}}`: Array of tuples containing deme (or nothing) and sampling time
"""
struct SampleConfiguration 
    sconf::Array{Tuple{Union{Nothing,String}, Float64}}
end 

"""
    SampleConfiguration(conffn::String)

Create a SampleConfiguration from a YAML configuration file.

# Arguments
- `conffn::String`: Path to the YAML configuration file

# Returns
- `SampleConfiguration`: The constructed sampling configuration
"""
function SampleConfiguration(conffn::String)
    SampleConfiguration(confstr = read(conffn, String))
end

"""
    SampleConfiguration(; confstr::String)

Create a SampleConfiguration from a YAML configuration string.

# Arguments
- `confstr::String`: YAML configuration string defining the sampling scheme

# Returns
- `SampleConfiguration`: The constructed sampling configuration
"""
function SampleConfiguration(; confstr::String)
    conf = YAML.load(confstr)
    samps = conf["sample"] 

    deme = Array{String}(undef,0)
    time = Array{Float64}(undef,0)

    for s in samps 
        if "table" ∈ keys(s)
            stx = DataFrame(CSV.File(s["table"]))
            push!(time, stx.sample_time...)
            demes = [String(x) for x in stx.deme] # may need to promote from String1
            push!(deme, demes...)
        elseif "time" ∈ keys(s)
            if s["time"] isa String
                tx = (Meta.parse(s["time"])) |> eval |> collect 
            elseif s["time"] isa Number
                tx = s["time"] 
            else 
                throw(s)
            end
            ltx = length(tx)
            if ltx == 1 && "size" ∈ keys(s)
                @assert s["size"] isa Number 
                ltx = s["size"] 
                tx = repeat([tx], ltx)
            end
            push!(time, tx...)
            if "deme" ∈ keys(s) 
                push!(deme, repeat([s["deme"]], ltx)...)
            else 
                push!(deme, repeat(["NULL"], ltx)...)
            end
        end
    end 
    SampleConfiguration(collect(zip(deme, time)))
end
