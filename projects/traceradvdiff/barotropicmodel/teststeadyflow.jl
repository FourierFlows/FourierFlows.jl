include("../../../src/fourierflows.jl")

import FourierFlows
import FourierFlows.TracerAdvDiff

n = 256
L = 2π
H = 1.0
κ = 1e-4

# Parameters
 ε = 0.2
k₁ = 2π/L      
 δ = L/n
 U = 1.0
um = U / (L*(1-ε))
tf = H*L/U

CFL = 0.2
Δt₁ = CFL*δ/um
nsteps₁ = tf/Δt₁
savefreq = 100

# Total steps, intermediate steps.
intsteps = ceil(Int, nsteps₁/savefreq) # Number of snapshots and diags
totsteps = intsteps*savefreq
Δt = tf/totsteps

# Sinsuidal topography and barotropic velocity
 h(x) = H(1 - ε*sin(k₁*x))       # Topography
∂h(x) = -H*k₁*ε*cos(k₁*x)        # Topographic x-gradient

u(x, y) = U/h(x)                 # Barotropic x-velocity
v(x, y) = y*u(x, y)*∂h(x)/h(x)   # Topography-following z-velocity

# Initial condition
x₀, y₀ = L/2, H/2
ci(x, y) = exp( -(x-x₀)^2/(2*δx^2) - (y-y₀)^2/(2*δy^2) ) / (2π*δx*δy)

# Initialize the problem
g = FourierFlows.TwoDGrid(n, L, n, H+ε; y0=-H-ε)
prob = TracerAdvDiff.ConstDiffSteadyFlowProblem(grid=g, κ=κ, u=u, v=v) 

# Output
getc(prob) = irfft(prob.vars.sol, prob.grid.nx)
outc = Output(prob, filename, (:c, getc))
outputs = [outc]

# Diagnostics
variance(prob) = Cy2(prob.vars.c, prob.grid)
variancediag = Diagnostic(variance, prob)
diags = [variancediag]


# Plotting
fig, axs = subplots(figsize=(8, 4))
mask = Array{Bool}(n, n)
for j=1:n, i=1:n;
  mask[i, j] = g.y[j] < -h(g.x[i]) ? true : false 
end

# Integrate
while prob.step < totsteps

  stepforward!(prob, diags; nsteps=substeps)
  pcolormesh(prob.grid.x, prob.grid.y, outputs[:c])
  pause(0.1)
  
end
