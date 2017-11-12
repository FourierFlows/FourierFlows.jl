include("./plottools.jl")

using PyPlot, PyCall

import SLSPlotTools, FourierFlows.TwoModeBoussinesq

import FourierFlows.TwoModeBoussinesq: mode0apv, mode1apv, mode1speed, mode1w

@pyimport mpl_toolkits.axes_grid1 as pltgrid

# Plot parameters
   simname = "nk32_nu1e+32_122f_n512_ep20_Ro05_nkw32"
plotprefix = "./plots/mitsls_"
  filename = "../data/$simname.jld2"

   eddylim = 20
       w0  = 0.08
     istep = 70




# Get problem and calculate σ
prob, steps = SLSPlotTools.load_boussinesq_ivp(filename)

nkw = 32
 kw = 2π*nkw/prob.grid.Lx
  α = (prob.params.N*kw)^2/(prob.params.f*prob.params.m)^2
  σ = prob.params.f*sqrt(1+α)
 tσ = 2π/σ
 Re = prob.grid.Lx/40




# Get snapshot
step = steps[istep]
x, y = prob.grid.X/Re, prob.grid.Y/Re
t, Zh, solc = SLSPlotTools.get_snapshot(filename, step)


# Analysis
Z = irfft(Zh, prob.grid.nx) 
u = ifft(solc[:, :, 1])
v = ifft(solc[:, :, 2])
p = ifft(solc[:, :, 3])

TwoModeBoussinesq.set_Z!(prob, Z)
TwoModeBoussinesq.set_uvp!(prob, u, v, p)

w = mode1w(prob)


# Plot
close("all")
fig, ax = subplots()

pcolormesh(x, y, w, cmap="RdBu_r", vmin=-w0, vmax=w0)

axis("equal")
ax[:set_adjustable]("box-forced")


ax[:set_xlim](-eddylim, eddylim)
ax[:set_ylim](-eddylim, eddylim)

ax[:axis]("off")
#ax[:tick_params](axis="both", which="both", length=0)
#ax[:xaxis][:set_ticks]([])
#ax[:yaxis][:set_ticks]([])

tight_layout()

savefig("refracted_wave.png", dpi=240, transparent=true)
