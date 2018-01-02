using Documenter, FourierFlows

makedocs(modules=[FourierFlows],
        doctest=true)

deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/FourierFlows/FourierFlows.jl.git",
    julia  = "0.6.0",
    osname = "linux")
