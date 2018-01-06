using Documenter, FourierFlows

makedocs(modules=[FourierFlows],
         format = :html,
         sitename = "navidccc",
         pages = Any[
          "Home" => "index.md",
          "Tutorials" => Any[
              "tutorials/page1.md",
              "tutorials/page2.md",
              "tutorials/page3.md"
          ],
          "section2" => Any[
              "sec2/page1.md",
              "sec2/page2.md",
              "sec2/page3.md"
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
