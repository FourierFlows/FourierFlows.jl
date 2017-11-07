include("./plottools.jl")

using PyPlot, PyCall

import SLSPlotTools, FourierFlows.TwoModeBoussinesq

import FourierFlows.TwoModeBoussinesq: mode0apv, mode1apv, mode1speed, mode1w, 
  wave_induced_speed, wave_induced_psi, wave_induced_uv, lagrangian_mean_uv,
  calc_chi, calc_chi_uv, totalenergy, mode0energy, mode1energy, CFL,
  lagrangian_mean_psi, apv_induced_psi

@pyimport mpl_toolkits.axes_grid1 as pltgrid

# Names
   simname = "nk32_nu1e+32_122f_n512_ep20_Ro05_nkw32"
plotprefix = "./plots/mitsls_"
  filename = "../data/$simname.jld2"

prob, steps = SLSPlotTools.load_boussinesq_ivp(filename)

# Plot parameters
eddylim = 20
    Ro0 = 0.04
    w0  = 0.08
    sp0 = 1.0
  istep = 80


# Non-retrievable parameters
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

Umax = maximum(sqrt.(prob.vars.U.^2+prob.vars.V.^2))

w = mode1w(prob)
Q = mode0apv(prob)
sp = wave_induced_speed(σ, prob)/Umax
Psiw = wave_induced_psi(σ, prob)
PsiL = lagrangian_mean_psi(σ, prob)
PsiQ = apv_induced_psi(prob)
Chi  = calc_chi(prob)
uL, vL = lagrangian_mean_uv(σ, prob)
uc, vc = calc_chi_uv(prob)

spc = sqrt.(uc.^2.0+vc.^2.0)/Umax




# Plot
close("all")
fig, ax = subplots()

pcolormesh(x, y, w, cmap="RdBu_r")

axis("equal")
ax[:set_adjustable]("box-forced")
ax[:set_xlim](-eddylim, eddylim)
ax[:set_ylim](-eddylim, eddylim)
ax[:tick_params](axis="both", which="both", length=0)

#ax[:set_xlabel](L"x/R")
#ax[:set_xlabel](L"x/R")
ax[:xaxis][:set_ticks]([])
ax[:yaxis][:set_ticks]([])

tight_layout()

savefig("refracted_wave.png", dpi=240)
