include("../src/fourierflows.jl")

using FourierFlows

import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy

# Physical parameters
 n = 128
 L = 2π
 nν = 4
  ν = 0e-8
withfilter = true

# Time-stepping
dt = 5e-3
nsteps = 8000
nsubs  = 200


# Initialize with random numbers
prob = TwoDTurb.InitialValueProblem(n, L, ν, nν, dt, withfilter)

# Initial condition closely following pyqg barotropic example
# that reproduces the results of the paper by McWilliams (1984)
srand(1234)
k0, E0 = 6, 0.5
modk = sqrt.(prob.grid.KKrsq)
psik = zeros(prob.grid.nk, prob.grid.nl)
psik =  (modk.^2 .* (1 + (modk/k0).^4)).^(-0.5)
psik[1, 1] = 0.0
psih = (randn(prob.grid.nkr, prob.grid.nl)+im*randn(prob.grid.nkr, prob.grid.nl)).*psik
psih = psih.*prob.ts.filter
Ein = real(sum(prob.grid.KKrsq.*abs2.(psih)/(prob.grid.nx*prob.grid.ny)^2))
psih = psih*sqrt(E0/Ein)
qi = -irfft(prob.grid.KKrsq.*psih, prob.grid.nx)
qh = rfft(qi)
TwoDTurb.set_q!(prob, qi)
TwoDTurb.updatevars!(prob)
prob.vars.psih = - prob.vars.qh .* prob.grid.invKKrsq

g = prob.grid
v = prob.vars

Q0_parseval2 = 0.5*FourierFlows.parsevalsum2(v.qh, g)
println("calculate enstrophy using parsevalsum2")
println("Q0 (parseval2) = ", Q0_parseval2)
println(" ")

Q0_parseval = 0.5*FourierFlows.parsevalsum(abs2.(v.qh), g)
println("calculate enstrophy using parsevalsum")
println("Q0 (parseval) = ", Q0_parseval)
println(" ")

Q0_integral = 0.5*sum(v.q.^2)*g.dx*g.dy
println("calculate enstrophy using integral")
println("Q0 (integral) = ", Q0_integral)
println(" ")

E0_parseval2 = 0.5*FourierFlows.parsevalsum2(sqrt.(g.KKrsq).*abs.(v.psih), g)
println("calculate enstrophy using parsevalsum2")
println("E0 (parseval2) = ", E0_parseval2)
println(" ")

E0_parseval = 0.5*FourierFlows.parsevalsum(g.KKrsq.*abs2.(v.psih), g)
println("calculate enstrophy using parsevalsum")
println("E0 (parseval) = ", E0_parseval)
println(" ")

E0_integral = 0.5*sum(v.U.^2 + v.V.^2)*g.dx*g.dy
println("calculate enstrophy using integral")
println("E0 (integral) = ", E0_integral)
println(" ")
