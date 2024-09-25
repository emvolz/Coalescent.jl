using Documenter
using Coalescent 
push!(LOAD_PATH,"../src/")

makedocs(
	sitename = "Coalescent.jl"
	# , modules = [Coalescent]
	# , format = Documenter.HTML()
	# , pages = [
	# 	"Home" => "README.jmd",  
	# ]
)


deploydocs( repo = "github.com/emvolz/Coalescent.jl.git" )
