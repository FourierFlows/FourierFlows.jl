include("rochadipole.jl")

import FourierFlows.NIWQG


# Setup
nx = 256
dtfrac = 0.0025
nsteps = 10000
substeps = 200

prob = makedipoleproblem(nx=nx, dtfrac=dtfrac,
  kape=5e7, nuw=1e7, Uw=1e-1)

function energyeqn(prob)
  NIWQG.updatevars!(prob)
  niwke      = NIWQG.niwke(prob)
  niwpe      = NIWQG.niwpe(prob)
  qgenergy   = NIWQG.qgenergy(prob)
  niwkediss  = NIWQG.niwke_dissipation(prob) 
  niwpediss  = NIWQG.niwpe_dissipation(prob) 
  qgdiss     = NIWQG.qg_dissipation(prob) 
  waveqgdiss = NIWQG.waveqg_dissipation(prob) 
  [niwke, niwpe, qgenergy, niwkediss, niwpediss, qgdiss, waveqgdiss]
end

ediags = Diagnostic(energyeqn, prob; nsteps=nsteps)


# Plotting
Ro = maximum(abs.(prob.vars.q)) / prob.params.f
sp = maximum(abs.(prob.vars.phi))

rossbyq(prob)   = prob.vars.q / prob.params.f
wavespeed(prob) = abs.(prob.vars.phi)
message(prob)   = @sprintf("\$t = %.1f\$ inertial periods", 
  prob.t*prob.params.f/(2*pi))


basicplot = FourierFlows.TwoComponentProblemPlot(prob,
  rossbyq,   [-Ro, Ro],   "RdBu_r",
  wavespeed, [0.0, 2*sp], "YlGnBu_r",
  message, "./plots/test")
  


function steppingmsg(prob, ediags) 

  keres = sum(ediags.value[1:2]) / ediags.data[1][1] 
  qgres = sum(ediags.value[3:7]) / sum(ediags.data[1][3:4]) 

  @sprintf(
    "step: %04d, niw ke: %.3f (res: %.2e), niw pe: %.3f, qg: %.3f (res: %.2e)",
    prob.step, ediags.value[1]/ediags.data[1][1], keres,
    ediags.value[2] / sum(ediags.data[1][2:3]),
    ediags.value[3] / sum(ediags.data[1][2:3]), qgres)
end

# Loop
@printf("Time-stepping. Ro: %.3f, Uw: %.3f", Ro, sp)
while prob.step < nsteps
  FourierFlows.stepforward!(prob, ediags; nsteps=substeps)
  NIWQG.updatevars!(prob)
  FourierFlows.makeplot!(basicplot; show=true, save=false)
  println(steppingmsg(prob, ediags))
end
