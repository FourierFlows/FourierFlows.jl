using Documenter, FourierFlows

makedocs(
   modules = [FourierFlows],
   doctest = false, clean = true,
    format = :html,
   authors = "Gregory L. Wagner and Navid C. Constantinou",
  sitename = "FourierFlows.jl",
     pages = Any[
              "Home" => "index.md",
              "Code Basics" => "basics.md",
              "Forcing" => "forcing.md",
              "Modules" => Any[
                "modules/kuramotosivashinsky.md",
                "modules/twodturb.md",
                "modules/barotropicqg.md",
                "modules/boussinesq.md",
                "modules/traceradvdiff.md"
              ],
              "DocStrings" => Any["man/docstrings.md"]
             ]
)

deploydocs(
       repo = "github.com/FourierFlows/FourierFlows.jl.git",
     target = "build",
      julia = "0.6.3",
     osname = "linux",
       deps = nothing,
       make = nothing
 )
