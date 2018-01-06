using Documenter, FourierFlows

makedocs(modules=[FourierFlows],
         format = :html,
         sitename = "FourierFlows.jl",
         pages = Any[
          "Home" => "index.md",
          "DocStrings" => Any[
              "man/docstrings.md"
              ]
          ])

deploydocs(
     repo   = "github.com/FourierFlows/FourierFlows.jl.git",
     target = "build",
     julia = "0.6.0",
     osname = "linux",
     deps   = nothing,
     make   = nothing
 )
