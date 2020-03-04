using
  Documenter,
  FourierFlows

makedocs(
    modules = [FourierFlows],
      clean = true,
    doctest = false,
  checkdocs = :all,
     format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Gregory L. Wagner and Navid C. Constantinou",
   sitename = "FourierFlows.jl",
  
      pages = Any[
              "Home" => "index.md",
              "Code Basics" => "basics.md",
              "Forcing" => "forcing.md",
              "DocStrings" => Any[
                  "man/types.md",
                  "man/functions.md"]
                 ]
)

deploydocs(
  repo = "github.com/FourierFlows/FourierFlows.jl.git",
)
