include("./setup_eddywave.jl")



# -- Parameters --
  nkw = 32    
    n = 512 
    L = 2π*100e3*nkw
    α = 0.5             # Frequency parameter
    ε = 2e-1            # Wave amplitude
   Ro = 5e-2            # Eddy Rossby number
Reddy = L/40           # Eddy radius

passiveapv = true


# Setup
ew = EddyWave(L, α, ε, Ro, Reddy; nkw=nkw, dtfrac=1e-2, nsubperiods=1,
  nν0=8, nν1=8, ν0=1e32, ν1=1e16) 

plotpath = "./plots"
plotname = @sprintf("nk32_nu%.0e", ew.ν0)
  
prob, diags = eddywavesetup(n, ew, perturbwavefield=false, 
  passiveapv=passiveapv)
etot, e0, e1 = diags[1], diags[2], diags[3]


# Output
fileprefix = @sprintf("./data/%s_%df_n%d_ep%02d_Ro%02d_nkw%02d",
  plotname, 100*ew.σ/ew.f, n, 100ε, 100Ro, ew.nkw)

i, testprefix = 0, fileprefix
while isfile(testprefix*".jld2"); i+=1; testprefix=fileprefix*"_$i"; end
filename = testprefix * ".jld2"

saveproblem(prob, filename)
saveeddywave(ew, filename)

getsolr(prob) = prob.vars.solr
getsolc(prob) = prob.vars.solc

outs = [
  Output("solr", getsolr, prob, filename),
  Output("solc", getsolc, prob, filename),
]


# Plotting
fig, axs = subplots(ncols=2, nrows=2, figsize=(8, 8)) 
xr, yr = prob.grid.X/Reddy, prob.grid.Y/Reddy


# Run
startwalltime = time()
while prob.step < ew.nsteps

  stepforward!(prob, diags; nsteps=ew.nsubs)
  TwoModeBoussinesq.updatevars!(prob)
  saveoutput(outs)

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

  plotmsg1 = @sprintf(
    "\$t=% 3d\$ wave periods, \$E_0=%.3f\$, ",
    round(Int, prob.t/ew.twave), e0.value/e0.data[1])
    
  plotmsg2 = @sprintf(
    "\$E_1=%.3f\$, \$E_{\\mathrm{tot}}=%.6f\$",
    e1.value/e1.data[1], etot.value/etot.data[1])

  savename = @sprintf("%s_n%d_%02df_ep%02d_Ro%02d_%06d.png", 
    joinpath(plotpath, plotname), n, floor(Int, 100ew.σ/ew.f), 
    floor(Int, 100*ew.ε), floor(Int, 100*ew.Ro), prob.step)

  close("all")
  fig, axs = subplots(ncols=2, nrows=2, figsize=(8, 8)) 
  makefourplot!(axs, prob, ew, xr, yr, savename; 
    message=plotmsg1*plotmsg2, save=true, eddylim=nothing, show=false,
    passiveapv=passiveapv)

end
