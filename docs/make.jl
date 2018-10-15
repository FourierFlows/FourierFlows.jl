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
              "Modules" => Any[
                "modules/kuramotosivashinsky.md",
                "modules/traceradvdiff.md",
              ],
              "DocStrings" => Any[
              "man/types.md",
              "man/functions.md"]
             ]
)

deploydocs(
       repo = "github.com/FourierFlows/FourierFlows.jl.git",
     target = "build",
      julia = "0.7",
     osname = "linux",
       deps = nothing,
       make = nothing
 )
