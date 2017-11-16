include("../../../src/fourierflows.jl")

n = 256
L = 2π
H = 1.0
κ = 1e-4

xi = 

u(x, y) = 1.0
v(x, y) = 0.0
ci(x, y) = exp( -(x-xi)^2/(2*δx^2) - (y-yi)^2/(2*δy^2) ) / (2π*δx*δy)

prob = ConstDiffSteadyFlowProblem(n, L, κ, u, v) 
