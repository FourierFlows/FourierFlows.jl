include("./setup_eddywave.jl")



# -- Parameters --
  nkw = 8
    n = 256
    L = 2π*100e3*nkw
    α = 0.5             # Frequency parameter
    ε = 1e-1            # Wave amplitude
   Ro = 2e-2            # Eddy Rossby number
Reddy = L/10            # Eddy radius


# Setup
ew = EddyWave(L, α, ε, Ro, Reddy; nkw=nkw, dtfrac=2e-2, nsubperiods=1,
  nν0=8, nν1=6, ν0=1e28, ν1=1e6) 

plotpath = "./plots"
plotname = @sprintf("nu0_%.0e", ew.ν0)
  
prob, diags = eddywavesetup(n, ew, perturbwavefield=false)
etot, e0, e1 = diags[1], diags[2], diags[3]


# Output
fileprefix = @sprintf("./data/%s_%df_n%d_ep%02d_Ro%02d_nkw%02d",
  plotname, 100*ew.σ/ew.f, n, 100ε, 100Ro, ew.nkw)

i, testprefix = 0, fileprefix
while isfile(testprefix*".jld2"); i+=1; testprefix=fileprefix*"_$i"; end
filename = testprefix * ".jld2"

saveproblem!(prob, filename)


getZ(prob)    = prob.vars.Z
getsolr(prob) = prob.vars.solr
getsolc(prob) = prob.vars.solc
getpsiw(prob) = wave_induced_psi(ew.σ, prob)
getuw(prob)   = wave_induced_uv(ew.σ, prob)[1]
getvw(prob)   = wave_induced_uv(ew.σ, prob)[2]


outs = [
  Output("Z",    getZ,     prob, filename),
  Output("solr", getsolr,  prob, filename),
  Output("solc", getsolc,  prob, filename),
  Output("Q",    mode0apv, prob, filename),
  Output("w",    mode1w,   prob, filename),
  Output("uw",   getuw,    prob, filename),
  Output("vw",   getvw,    prob, filename)
]



# Plotting
fig, axs = subplots(ncols=3, nrows=1, figsize=(12, 5)) 
xr, yr = prob.grid.X/Reddy, prob.grid.Y/Reddy

tperiods = 0

savename = @sprintf("%s_n%d_%02df_ep%02d_Ro%02d_%06d.png", 
  joinpath(plotpath, plotname), n, floor(Int, 100ew.σ/ew.f), 
  floor(Int, 100*ew.ε), floor(Int, 100*ew.Ro), prob.step)

cbs, msg = makethreeplot!(axs, prob, ew, xr, yr, savename; 
  message="$tperiods", save=true, eddylim=nothing, show=false, 
  makecolorbar=true)

# Run
startwalltime = time()
while prob.step < ew.nsteps

  stepforward!(prob, diags; nsteps=ew.nsubs)
  TwoModeBoussinesq.updatevars!(prob)
  saveoutput!(outs)

  walltime = (time()-startwalltime)/60

  log1 = @sprintf(
    "step: %04d, t: %d, CFL: %.2f, max Ro: %.4f, ", 
    prob.step, prob.t/ew.twave, CFL(prob, prob.ts.dt), 
    maximum(abs.(prob.vars.Z))/ew.f)

  log2 = @sprintf(
    "ΔE: %.3f, Δe: %.3f, Δ(E+e): %.6f, τ: %.2f min",
    e0.value/e0.data[1], e1.value/e1.data[1], etot.value/etot.data[1],
    walltime)

  println(log1*log2)

  savename = @sprintf("%s_n%d_%02df_ep%02d_Ro%02d_%06d.png", 
    joinpath(plotpath, plotname), n, floor(Int, 100ew.σ/ew.f), 
    floor(Int, 100*ew.ε), floor(Int, 100*ew.Ro), prob.step)

  tperiods = round(Int, prob.t/ew.twave)
  cbs, msg = makethreeplot!(axs, prob, ew, xr, yr, savename; 
    message="$tperiods", save=true, eddylim=nothing, show=false, 
    removemsg=msg)


end
