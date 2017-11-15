include("../../../src/fourierflows.jl")

using FourierFlows, PyPlot, PyCall, JLD2

import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy

@pyimport mpl_toolkits.axes_grid1 as pltgrid

L  = 2π*1600e3
n  = 512
nν = 4
ν  = 1e12
f  = 1e-4
α  = 1.0 #(N*kw)^2/(f*m)^2
σ  = f*sqrt(1+α)
tσ = 2π/σ
dt = 5e-2 * tσ

nperiods    = 400
nsubperiods = 1
nsteps = ceil(Int, nperiods*tσ/dt)
nsubs  = ceil(Int, nsubperiods*tσ/dt)

savename = "./data/twodturb_n512_Ro20_nnu4_nu1e+12.jld2"

@load savename Z

prob = TwoDTurb.InitialValueProblem(n, L, ν, nν, dt)
TwoDTurb.set_q!(prob, Z)

# Diagnostics
E = Diagnostic(energy, prob; nsteps=nsteps)

# Prepare output
fileprefix = savename[1:end-5] * "_fwd"
plotprefix = savename[8:end-5] * "_fwd"

i, testprefix = 1, fileprefix
while isfile(testprefix*".jld2"); i+=1; testprefix=fileprefix*"-$i"; end
filename = testprefix * ".jld2"

getsol(prob) = prob.vars.sol
output = Output(prob, filename, (:sol, getsol))

startwalltime, iplot = time(), 0
while prob.step < nsteps

  stepforward!(prob, E; nsteps=nsubs)
  saveoutput(output)

  log = @sprintf("step: %04d, t: %d, ΔE: %.6f, τ: %.2f min", 
    prob.step, prob.t/tσ, E.value/E.data[1], (time()-startwalltime)/60)
    
  println(log)

  if prob.step/nsubs % 10 == 0
    TwoDTurb.updatevars!(prob) 
    close("all")
    fig, axs = subplots()
    pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
    axis("equal")
    title(@sprintf("max q/f: %.2f", maximum(abs.(prob.vars.q)/f)))

    plotname = @sprintf("./plots/%s_%04d.png", plotprefix, iplot)
    savefig(plotname, dpi=240)
    iplot += 1
  end

end
