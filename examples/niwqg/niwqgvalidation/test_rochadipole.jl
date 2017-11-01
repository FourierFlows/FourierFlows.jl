include("rochadipole.jl")

#import FourierFlows


# Setup
nx = 256
dtfrac = 0.05
nsteps = 10000
substeps = round(Int, 1.0/dtfrac)

prob = make_dipole_problem(nx=nx, dtfrac=dtfrac,
  kape=1e8, nuw=3e8, Uw=2e-1)




# Plotting
Ro = maximum(abs.(prob.vars.q)) / prob.params.f
sp = maximum(abs.(prob.vars.phi))

rossbyq(prob)   = prob.vars.q / prob.params.f
wavespeed(prob) = abs.(prob.vars.phi)
message(prob)   = @sprintf("\$t = %.1f\$ inertial periods", 
  prob.t*prob.params.f/(2*pi))

basicplot = FourierFlows.TwoComponentProblemPlot(prob,
  rossbyq,   [-Ro, Ro], "RdBu_r",
  wavespeed, [0.0, 2*sp], "YlGnBu_r",
  message, "./plots/test")
  



# Loop
@printf("Time-stepping. Ro: %.3f, Uw: %.3f", Ro, sp)
while prob.step < nsteps
  FourierFlows.stepforward!(prob; nsteps=substeps)
  NIWQG.updatevars!(prob)

  FourierFlows.makeplot!(basicplot; show=false, save=true)

  println(prob.step)
end
