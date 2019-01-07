using
  Documenter,
  FourierFlows

makedocs(
   modules = [FourierFlows],
   doctest = false, clean = true,
 checkdocs = :all,
    format = :html,
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
     target = "build",
       deps = nothing,
       make = nothing
 )
