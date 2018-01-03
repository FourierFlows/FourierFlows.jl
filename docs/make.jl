using Documenter, FourierFlows

makedocs(modules=[FourierFlows],
        doctest=true,
        format = :html,
        sitename = "FourierFlows.jl",
        authors="Gregory L. Wagner and Navid C. Constantinou",
        pages = [
                 "Home" => "index.md"
                 ])

deploydocs(
    repo = "github.com/FourierFlows/FourierFlows.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)

# deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
#     repo = "github.com/FourierFlows/FourierFlows.jl.git",
#     julia  = "0.6.0",
#     osname = "linux")
