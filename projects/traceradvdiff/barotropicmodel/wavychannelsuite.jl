include("./periodicstraintools.jl")

name = "try_1"
 κ = 1e-4
 η = 1e-5
δx = 0.02
δy = 0.04


for ε in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
  wavybarotropicrun(ε, name=name, κ=κ, η=η, δx=δx, δy=δy)
end
