include("./periodicstraintools.jl")

name = "try_1"
 n =  256
 κ = 1e-4
 η = 1e-4
δx = 0.01
δy = 0.01
CFL = 0.01

ε = 0.2
wavybarotropicrun(ε, n=n, name=name, κ=κ, η=η, δx=δx, δy=δy, CFL=CFL)

#for ε in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
#  wavybarotropicrun(ε, n=n, name=name, κ=κ, η=η, δx=δx, δy=δy, CFL=CFL)
#end
