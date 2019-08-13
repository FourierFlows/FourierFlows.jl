# Workaround for JuliaLang/julia/pull/28625
if Base.HOME_PROJECT[] !== nothing
  Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using
  Documenter,
  FourierFlows

makedocs(
   modules = [FourierFlows],
   doctest = false,
     clean = true,
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

deploydocs(repo = "github.com/FourierFlows/FourierFlows.jl.git")
