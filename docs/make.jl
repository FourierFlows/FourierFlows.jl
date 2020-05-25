using
  Documenter,
  FourierFlows

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
      canonical = "https://fourierflows.github.io/FourierFlowsDocumentations.jl/dev/"
)

makedocs(
    modules = [FourierFlows],
    doctest = false,
      clean = true,
   checkdocs = :all,
     format = format,
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

deploydocs(        repo = "github.com/FourierFlows/FourierFlowsDocumentations.jl.git",
               versions = ["stable" => "v^", "v#.#"],
           push_preview = true,
)
