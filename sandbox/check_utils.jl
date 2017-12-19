include("../src/fourierflows.jl")

using FourierFlows

import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy

# Physical parameters
 nx = 128
 ny = 256
 Lx = 2π
 Ly = 3π
 nν = 4
  ν = 0e-8
withfilter = true

# Time-stepping
dt = 5e-3
nsteps = 8000
nsubs  = 200


# Initialize with random numbers
prob = TwoDTurb.InitialValueProblem(nx, Lx, ny, Ly, ν, nν, dt, withfilter)

g = prob.grid
v = prob.vars

# Simpler initial condition
# with E0 = 0.5*(m^2+n^2) and Q0 =0.5*(m^2+n^2)^2
m, n = 1*(2π/Lx), 2*(2π/Ly)
qi = -(m^2+n^2)*sin.(m*g.X).*sin.(n*g.Y)

TwoDTurb.set_q!(prob, qi)
TwoDTurb.updatevars!(prob)
v.psih = - v.qh .* g.invKKrsq

Q0_parseval2 = 0.5*FourierFlows.parsevalsum2(v.qh, g) / (Lx*Ly)
println("calculate enstrophy using parsevalsum2")
println("Q0 (parseval2) = ", Q0_parseval2)
println(" ")

Q0_parseval = 0.5*FourierFlows.parsevalsum(abs2.(v.qh), g) / (Lx*Ly)
println("calculate enstrophy using parsevalsum")
println("Q0 (parseval) = ", Q0_parseval)
println(" ")

Q0_integral = 0.5*sum(v.q.^2)*g.dx*g.dy / (Lx*Ly)
println("calculate enstrophy using integral")
println("Q0 (integral) = ", Q0_integral)
println(" ")

E0_parseval2 = 0.5*FourierFlows.parsevalsum2(sqrt.(g.KKrsq).*abs.(v.psih), g) / (Lx*Ly)
println("calculate enstrophy using parsevalsum2")
println("E0 (parseval2) = ", E0_parseval2)
println(" ")

E0_parseval = 0.5*FourierFlows.parsevalsum(g.KKrsq.*abs2.(v.psih), g) / (Lx*Ly)
println("calculate enstrophy using parsevalsum")
println("E0 (parseval) = ", E0_parseval)
println(" ")

E0_integral = 0.5*sum(v.U.^2 + v.V.^2)*g.dx*g.dy / (Lx*Ly)
println("calculate enstrophy using integral")
println("E0 (integral) = ", E0_integral)
println(" ")
