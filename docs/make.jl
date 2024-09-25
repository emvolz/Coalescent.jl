using Documenter
using Coalescent 
push!(LOAD_PATH,"../src/")

makedocs(;
	sitename = "Coalescent.jl"
	, modules = [Coalescent]
	, format = Documenter.HTML()
	, pages = [
		"Introduction" => "intro.md" 
	 	, "Reference" => "index.md",  
	]
)


deploydocs( repo = "github.com/emvolz/Coalescent.jl.git" )
